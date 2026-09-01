# CLAUDE.md

Repository-specific guidance for Claude Code sessions working in this repository.

## Provenance rule

Any factual claim in this repository's prose (README, CITATION.cff,
commit messages, docstrings) about scope, application, provenance,
dates, or instrument context MUST be traceable to a specific file
and line within THIS repository — never inferred from another
project worked on in the same session, even one involving the same
author or a related instrument. If no such source exists in this
repository, mark the claim as "[NEEDS CONFIRMATION FROM USER]"
instead of writing a plausible-sounding substitute. Before adding
any such claim, state which file/line justifies it.

## Why this rule exists

This repository's first published version incorrectly attributed its
application context to the TCABR tokamak, carried over from a
separate, unrelated project (tcabr-photon-budget) worked on in the
same session, involving the same author. The actual application,
stated in the wavelength-calibration scripts' own headers, is the
spectral characterization of a cold dielectric-barrier-discharge
(DBD) plasma. The error was caught after publication and required a
history rewrite (squash to a single corrected commit, tag v1.0.0
force-pushed) to remove it from the public record. This rule exists
to prevent a repeat.
