{% extends 'markdown/index.md.j2'%}

{% block input %}
{%- set text = cell.source| replace("\nshow()", "") %}
{%- if text %}
```julia
{%- if text.endswith(";") %}
{{ text[:-1] }}
{%- else %}
{{ text }}
{%- endif %}
```
{% endif %}
{% endblock input %}

{%- block markdowncell %}
{#- Jupyter/MathJax renders bare LaTeX display environments (e.g. a top-level
    \begin{equation}...\end{equation}) as math, but Documenter/KaTeX only picks
    up math delimited by $, $$ or ```math fences. Normalise those bare display
    environments to $$...$$ so they flow through the handling below, otherwise
    their contents leak into the prose and the underscores get mangled as
    Markdown emphasis. Note: environments meant to live *inside* math (split,
    array, aligned, cases, ...) are intentionally left untouched. -#}
{%- set src = cell.source
      | replace("\\begin{equation*}", "$$") | replace("\\end{equation*}", "$$")
      | replace("\\begin{equation}", "$$")  | replace("\\end{equation}", "$$")  %}
{%- set parts = src.split("$$") %}
{%- for i in range(parts | count) %}
{%- if i % 2 %}
```math
{{ parts[i] }}
```
{%- else %}{{ parts[i] | replace("$", "``") }}{%- endif %}
{%- endfor %}
{% endblock markdowncell %}
