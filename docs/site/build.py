#!/usr/bin/env python3
"""Build the static documentation index for tpt-dynamics.

Converts the project's Markdown documentation into HTML pages under docs/site/
(the same directory as this script) and writes docs/site/index.html and
docs/site/style.css. The output is a local, static documentation index. It is
not a hosted site and nothing is published by this script.

Standard library only. Run from anywhere:

    python docs/site/build.py

The converter handles the Markdown used in this repository: ATX headings,
paragraphs, ordered and unordered lists (one level, nested items are
flattened), fenced code, pipe tables, block quotes, emphasis, code spans and
links. It is not a full CommonMark implementation.
"""

from __future__ import annotations

import html
import posixpath
import re
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
REPO_URL = "https://github.com/tpt-solutions/tpt-dynamics"

# (section title, directory relative to the repo root, glob)
SECTIONS = [
    ("Tutorials", "docs/tutorials", "*.md"),
    ("Architecture decision records", "docs/adr", "*.md"),
    ("Benchmarks", "docs/benchmarks", "*.md"),
    ("Migration guides", "docs/migration", "*.md"),
    ("Audit reports", "docs/audit", "*.md"),
]
# Single Markdown files directly under docs/.
TOP_LEVEL = ["docs/dependencies.md", "docs/fuzzing.md", "docs/release-checklist.md"]
# Rustdoc crate indexes. They exist only after `cargo doc --workspace --no-deps`.
RUSTDOC_CRATES = [
    "tpt_dynamics_core",
    "tpt_dynamics_data",
    "tpt_dynamics_compute",
    "tpt_dynamics_sim",
    "tpt_dynamics_cli",
]

TOKEN_RE = re.compile(r"`([^`]+)`|\[([^\]]+)\]\(([^)\s]+)\)")
CODE_RE = re.compile(r"`([^`]+)`")
HEADING_RE = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")
RULE_RE = re.compile(r"^(?:-{3,}|\*{3,}|_{3,})$")
LIST_RE = re.compile(r"^(?:[-*+]|\d+\.)\s+")
TABLE_SEP_RE = re.compile(r"^\|?\s*:?-+:?\s*(\|\s*:?-+:?\s*)*\|?$")
SCHEME_RE = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.\-]*:")


@dataclass
class Page:
    src: str  # repo-relative POSIX path, e.g. docs/adr/ADR-005-licence-policy.md
    out: str  # path relative to docs/site/, e.g. adr/ADR-005-licence-policy.html
    section: str
    title: str = ""
    body: str = ""


@dataclass
class Ctx:
    src: str
    out: str
    pages: dict[str, str] = field(default_factory=dict)


def slugify(text: str) -> str:
    t = re.sub(r"[`*_]", "", text).lower()
    t = re.sub(r"[^a-z0-9\s-]", "", t)
    return re.sub(r"\s+", "-", t.strip())


def emphasis(s: str) -> str:
    s = re.sub(r"\*\*(?=\S)(.+?)(?<=\S)\*\*", r"<strong>\1</strong>", s)
    return re.sub(r"(?<![\w*])\*(?=\S)(.+?)(?<=\S)\*(?![\w*])", r"<em>\1</em>", s)


def fmt(text: str) -> str:
    """Code spans and emphasis only (no links)."""
    parts = []
    pos = 0
    for m in CODE_RE.finditer(text):
        parts.append(emphasis(html.escape(text[pos : m.start()])))
        parts.append("<code>" + html.escape(m.group(1)) + "</code>")
        pos = m.end()
    parts.append(emphasis(html.escape(text[pos:])))
    return "".join(parts)


def rewrite_href(url: str, ctx: Ctx) -> str:
    if url.startswith("#") or SCHEME_RE.match(url):
        return url
    path, _, anchor = url.partition("#")
    anchor = "#" + anchor if anchor else ""
    resolved = posixpath.normpath(posixpath.join(posixpath.dirname(ctx.src), path))
    if path.endswith(".md") and resolved in ctx.pages:
        target = ctx.pages[resolved]
        base = posixpath.dirname(ctx.out) or "."
        return posixpath.relpath(target, base) + anchor
    # Anything else (source files, directories) points at the repository on GitHub.
    kind = "tree" if (REPO / resolved).is_dir() else "blob"
    return f"{REPO_URL}/{kind}/master/{resolved}{anchor}"


def inline(text: str, ctx: Ctx) -> str:
    out = []
    pos = 0
    for m in TOKEN_RE.finditer(text):
        out.append(fmt(text[pos : m.start()]))
        if m.group(1) is not None:
            out.append("<code>" + html.escape(m.group(1)) + "</code>")
        else:
            href = html.escape(rewrite_href(m.group(3), ctx), quote=True)
            out.append(f'<a href="{href}">{fmt(m.group(2))}</a>')
        pos = m.end()
    out.append(fmt(text[pos:]))
    return "".join(out)


def split_row(line: str) -> list[str]:
    s = line.strip()
    if s.startswith("|"):
        s = s[1:]
    if s.endswith("|"):
        s = s[:-1]
    return [c.strip() for c in s.split("|")]


def is_block_start(s: str) -> bool:
    return bool(
        s.startswith("```")
        or HEADING_RE.match(s)
        or RULE_RE.match(s)
        or s.startswith(">")
        or LIST_RE.match(s)
        or s.startswith("|")
    )


def render_blocks(lines: list[str], ctx: Ctx) -> str:
    out: list[str] = []
    i = 0
    n = len(lines)
    while i < n:
        raw = lines[i]
        s = raw.strip()
        if not s:
            i += 1
            continue

        if s.startswith("```"):
            lang = s[3:].strip()
            j = i + 1
            code = []
            while j < n and not lines[j].strip().startswith("```"):
                code.append(lines[j])
                j += 1
            cls = f' class="language-{html.escape(lang, quote=True)}"' if lang else ""
            out.append(f"<pre><code{cls}>{html.escape(chr(10).join(code))}</code></pre>")
            i = j + 1
            continue

        m = HEADING_RE.match(s)
        if m:
            level = len(m.group(1))
            text = m.group(2)
            out.append(f'<h{level} id="{slugify(text)}">{inline(text, ctx)}</h{level}>')
            i += 1
            continue

        if RULE_RE.match(s):
            out.append("<hr>")
            i += 1
            continue

        if s.startswith("|") and i + 1 < n and TABLE_SEP_RE.match(lines[i + 1].strip()):
            head = split_row(s)
            i += 2
            rows = []
            while i < n and lines[i].strip().startswith("|"):
                rows.append(split_row(lines[i]))
                i += 1
            t = ['<div class="table-wrap"><table>', "<thead><tr>"]
            t += [f"<th>{inline(c, ctx)}</th>" for c in head]
            t += ["</tr></thead>", "<tbody>"]
            for r in rows:
                t.append("<tr>" + "".join(f"<td>{inline(c, ctx)}</td>" for c in r) + "</tr>")
            t.append("</tbody></table></div>")
            out.append("".join(t))
            continue

        if LIST_RE.match(s):
            ordered = s[0].isdigit()
            items: list[str] = []
            while i < n:
                cur = lines[i]
                cs = cur.strip()
                if not cs:
                    k = i + 1
                    while k < n and not lines[k].strip():
                        k += 1
                    if k < n and (lines[k][:1].isspace() or LIST_RE.match(lines[k].strip())):
                        i = k
                        continue
                    break
                indented = cur[:1].isspace()
                if LIST_RE.match(cs) and not indented:
                    items.append(LIST_RE.sub("", cs, count=1))
                    i += 1
                elif indented and items:
                    # Continuation text, or a nested item (flattened into the parent).
                    items[-1] += " " + cs
                    i += 1
                else:
                    break
            tag = "ol" if ordered else "ul"
            out.append(f"<{tag}>" + "".join(f"<li>{inline(t, ctx)}</li>" for t in items) + f"</{tag}>")
            continue

        if s.startswith(">"):
            quoted = []
            while i < n and lines[i].strip().startswith(">"):
                quoted.append(lines[i].strip()[1:].lstrip())
                i += 1
            out.append("<blockquote>" + render_blocks(quoted, ctx) + "</blockquote>")
            continue

        para = []
        while i < n and lines[i].strip() and not is_block_start(lines[i].strip()):
            para.append(lines[i].strip())
            i += 1
        if not para:  # unhandled block-start line: emit it as text so the loop advances
            para = [s]
            i += 1
        out.append("<p>" + inline(" ".join(para), ctx) + "</p>")
    return "\n".join(out)


def first_heading(text: str, fallback: str) -> str:
    for line in text.splitlines():
        m = re.match(r"^#\s+(.*)$", line.strip())
        if m:
            return re.sub(r"[`*_]", "", m.group(1)).strip()
    return fallback


def rel(target: str, from_out: str) -> str:
    return posixpath.relpath(target, posixpath.dirname(from_out) or ".")


def page_html(title: str, body: str, out: str, nav: list[tuple[str, list[Page]]], footer: str) -> str:
    nav_html = []
    for section, pages in nav:
        items = "".join(
            f'<li><a href="{rel(p.out, out)}">{html.escape(p.title)}</a></li>' for p in pages
        )
        nav_html.append(f"<h2>{html.escape(section)}</h2><ul>{items}</ul>")
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)} - tpt-dynamics docs</title>
<link rel="stylesheet" href="{rel("style.css", out)}">
</head>
<body>
<header><a href="{rel("index.html", out)}">tpt-dynamics documentation index</a></header>
<div class="layout">
<nav aria-label="Documentation">
<ul><li><a href="{rel("index.html", out)}">Index</a></li></ul>
{"".join(nav_html)}
</nav>
<main>
{body}
<footer>{footer}</footer>
</main>
</div>
</body>
</html>
"""


STYLE_CSS = """/* Generated by docs/site/build.py. Colours are tokens; dark mode follows the
   viewer's system setting unless data-theme is set on <html>. */
:root {
  --bg: #ffffff; --fg: #1d2330; --muted: #5b6474; --accent: #1f5fbf;
  --code-bg: #f3f5f8; --border: #d9dee7;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #0f1218; --fg: #e6e9ef; --muted: #9aa4b2; --accent: #7fb2ff;
    --code-bg: #1a1f2a; --border: #2a3140;
  }
}
:root[data-theme="dark"] {
  --bg: #0f1218; --fg: #e6e9ef; --muted: #9aa4b2; --accent: #7fb2ff;
  --code-bg: #1a1f2a; --border: #2a3140;
}
body { margin: 0; background: var(--bg); color: var(--fg);
  font: 16px/1.55 system-ui, -apple-system, "Segoe UI", sans-serif; }
header { border-bottom: 1px solid var(--border); padding: 12px 16px; }
header a { color: var(--fg); font-weight: 600; text-decoration: none; }
.layout { display: grid; grid-template-columns: 260px minmax(0, 1fr); gap: 24px;
  max-width: 1200px; margin: 0 auto; padding: 16px; }
nav ul { list-style: none; padding: 0; margin: 0; }
nav li { margin: 2px 0; }
nav h2 { font-size: 0.8rem; text-transform: uppercase; color: var(--muted);
  margin: 16px 0 4px; }
main { min-width: 0; max-width: 860px; }
a { color: var(--accent); }
pre { background: var(--code-bg); padding: 12px; overflow-x: auto; border-radius: 6px; }
code { font-family: ui-monospace, Consolas, "Courier New", monospace; font-size: 0.92em; }
:not(pre) > code { background: var(--code-bg); padding: 0 4px; border-radius: 4px; }
.table-wrap { overflow-x: auto; }
table { border-collapse: collapse; }
th, td { border: 1px solid var(--border); padding: 6px 8px; text-align: left; vertical-align: top; }
blockquote { margin: 0; padding: 0 12px; border-left: 3px solid var(--border); color: var(--muted); }
.note { border: 1px solid var(--border); border-radius: 6px; padding: 10px 12px; color: var(--muted); }
footer { margin-top: 32px; font-size: 0.85rem; color: var(--muted);
  border-top: 1px solid var(--border); padding-top: 8px; }
@media (max-width: 720px) {
  .layout { grid-template-columns: 1fr; padding: 12px 16px; }
}
"""


def collect_pages() -> list[Page]:
    pages: list[Page] = []
    for section, rel_dir, pattern in SECTIONS:
        d = REPO / rel_dir
        if not d.is_dir():
            continue
        for p in sorted(d.glob(pattern)):
            src = p.relative_to(REPO).as_posix()
            pages.append(Page(src, src[len("docs/") : -3] + ".html", section))
    for src in TOP_LEVEL:
        if (REPO / src).is_file():
            pages.append(Page(src, src[len("docs/") : -3] + ".html", "Project documentation"))
    return pages


def main() -> None:
    pages = collect_pages()
    page_map = {p.src: p.out for p in pages}

    # Render every Markdown page first so the navigation can show real titles.
    for p in pages:
        text = (REPO / p.src).read_text(encoding="utf-8")
        p.title = first_heading(text, Path(p.src).stem)
        ctx = Ctx(src=p.src, out=p.out, pages=page_map)
        p.body = render_blocks(text.splitlines(), ctx)

    nav: list[tuple[str, list[Page]]] = []
    for section in dict.fromkeys(p.section for p in pages):
        nav.append((section, [p for p in pages if p.section == section]))

    # Remove previously generated pages so renamed sources do not leave stale HTML.
    for old in HERE.rglob("*.html"):
        old.unlink()

    footer_tpl = 'Source: <a href="{url}">{src}</a> (rendered by docs/site/build.py)'
    for p in pages:
        footer = footer_tpl.format(url=f"{REPO_URL}/blob/master/{p.src}", src=p.src)
        target = HERE / p.out
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(page_html(p.title, p.body, p.out, nav, footer), encoding="utf-8")

    index_body = [
        "<h1>tpt-dynamics documentation index</h1>",
        '<p class="note">This is a <strong>static documentation index</strong>, generated '
        "from the repository's Markdown by <code>docs/site/build.py</code>. It is not a hosted "
        "documentation site. Nothing here is published. The links work when this folder is "
        "opened from a local checkout, or the same files can be read on GitHub. The Markdown "
        "files remain the source of truth.</p>",
    ]
    for section in dict.fromkeys(p.section for p in pages):
        items = "".join(
            f'<li><a href="{rel(p.out, "index.html")}">{html.escape(p.title)}</a>'
            f' <span class="note">({html.escape(p.src)})</span></li>'
            for p in pages
            if p.section == section
        )
        index_body.append(f"<h2>{html.escape(section)}</h2><ul>{items}</ul>")

    rustdoc_items = "".join(
        f'<li><a href="../../target/doc/{c}/index.html">{c}</a></li>' for c in RUSTDOC_CRATES
    )
    index_body.append(
        "<h2>API reference (rustdoc)</h2>"
        '<p class="note">Generated by <code>cargo doc --workspace --no-deps</code>, which CI also '
        "runs. These links point at the local <code>target/doc</code> output and are broken "
        "until that command has been run from the repository root.</p>"
        f"<ul>{rustdoc_items}</ul>"
    )
    index_body.append(
        "<h2>Regenerate</h2>"
        '<p class="note">From the repository root: <code>python docs/site/build.py</code>. '
        "Standard library only.</p>"
    )
    footer = "Static documentation index. Not a hosted site; nothing is published from this folder."
    (HERE / "index.html").write_text(
        page_html("Index", "\n".join(index_body), "index.html", nav, footer), encoding="utf-8"
    )
    (HERE / "style.css").write_text(STYLE_CSS, encoding="utf-8")
    print(f"wrote {len(pages)} pages and index.html under {HERE}")


if __name__ == "__main__":
    main()
