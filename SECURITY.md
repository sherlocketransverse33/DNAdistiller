# Security policy

## Reporting

Report vulnerabilities through [GitHub's private advisory form](https://github.com/Lucaks3/DNAdistiller/security/advisories/new). Please do not open a public issue for anything that would expose users before it is fixed.

**Never include a real DNA file in a report.** If reproducing a bug needs data, write a synthetic file with invented genotypes, or send the first ten lines with the rsIDs removed. A bug report is not worth your genome, and we cannot un-see a file once it is in an issue thread.

## What counts as a vulnerability here

This tool reads the most sensitive file most people own, so the threat model is narrower and stranger than for most software.

Treat as a security bug:

- **Any network call.** The core claim is that data cannot leave the machine because no code can send it. A socket opened anywhere in `src/dnadistiller/`, including in a dependency we pulled in, breaks that claim for every user. Report it.
- **Genetic data written to disk** anywhere the user did not name with `--out`. Temp files, caches, logs, and crash dumps all count. The original prototype this project replaced wrote uploads to a fixed path and leaked them on error, which is exactly the shape of the problem.
- **Redaction failures.** An rsID or genotype appearing in output at a tier that does not permit it. `tests/test_profile.py` guards this, and a case it misses is a real finding.
- **Genotypes in stack traces or error messages.** An exception that prints a user's call into a terminal, a log, or a crash report has disclosed it.
- **Catalog code execution.** Marker files are data. They are loaded with `yaml.safe_load` specifically so a YAML file cannot construct Python objects. A path that reaches `yaml.load`, `eval`, or `pickle` is a vulnerability.

Not security bugs, though still worth reporting as issues:

- A marker with the wrong rsID, effect allele, or strand. Serious and worth fixing quickly, but it is a correctness bug, not a disclosure one.
- Parsing failures on a provider's file.

## Design constraints

The project holds these regardless of what a feature request wants:

1. **No network code.** There is no analytics, no update check, no telemetry, and no remote catalog fetch. `pyproject.toml` has three runtime dependencies so that a user can plausibly audit the whole stack.
2. **The full genome never touches disk.** `Genome` lives in memory and is never serialised. Only `Profile`, which is redacted by construction, can be written out.
3. **Redaction happens in one function.** Every output format goes through `profile.redact`. A renderer cannot disclose a field it was never handed. Adding a format that reads `Finding` directly defeats this.
4. **Opt in, never opt out.** Filters build up from nothing. A bug in an opt-out filter discloses more than the user asked for; the same bug in an opt-in filter discloses less.

## Scope

The tool runs locally under your own account. It does not defend against someone who already controls your machine: if an attacker can read your home directory, they can read the DNA file directly and do not need us.

What it defends against is the ordinary failure of sending too much to someone else, on purpose, without realising how much it was.

## Supported versions

Pre-1.0. Fixes land on `main` only.
