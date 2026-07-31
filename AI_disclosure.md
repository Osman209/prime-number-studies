# AI assistance disclosure 

The wording below is deliberately specific about the division of labour, because a
vague acknowledgement is worse than none: a reader needs to know which parts to
check harder, not merely that a tool was involved.

---

## For the README

Put it under its own heading, near the reproducibility section.

### On the use of AI assistance

The verification scripts in `code/` and `harness/`, and much of the prose in
`papers/`, were written with the assistance of **Claude (Anthropic)**, used as a
working collaborator throughout: drafting and rewriting code, running the
computations, drafting and editing text, searching the literature, and — most
usefully — auditing the papers against their own scripts.

The research direction, the questions asked, the decisions about what to publish
and what to withdraw, and the final responsibility for every claim are the
author's.

A great many of the corrections recorded in these papers were found by that
auditing: a claim stated more precisely than the computation supported, a script
that had not caught up with a correction to the text, a formula quoted but never
tested. Several were errors the assistant had itself introduced and then found on
a later pass. Where a result is reported here, it is because a script regenerates
it and the script has been read; that discipline, rather than any assurance about
the tool, is what the reader is asked to rely on.

---

## For each paper

One line, placed immediately before the references.

**AI assistance.** The verification script accompanying this paper, and much of
its prose, were written with the assistance of Claude (Anthropic). The research
direction, the decisions, and the responsibility for every claim are the author's.
See the repository README for a fuller statement.

---

## Notes on placement

- Put the paper line **before** the references, not in a footnote — the point is
  that a reader sees it while deciding how much to trust the tables.
- The four existing papers and the two notes should all carry it, not just the new
  ones. Retrofitting it is a correction like any other and belongs in the same
  release.
- If a journal or preprint server asks for a disclosure in a particular form, use
  theirs; this is written for the repository.
