# Test fixtures

**Every file in this directory is synthetic. Never add real genetic data — not
yours, not a relative's, not a friend's, not one found in a public dataset.**

## Why this rule has no exceptions

Genetic data cannot be revoked. A password can be rotated and a leaked document
loses relevance, but a genome is a permanent identifier that stays valid for a
lifetime. Git makes that worse: a DNA file committed and then deleted in a later
commit is still in the history, still in every clone, and still on GitHub's
servers.

The part people miss is that it is not only your call to make. Your genome
discloses a great deal about your parents, siblings, and children — people who
never agreed to appear in a public repository, and who cannot opt out after the
fact.

"Data from a public dataset" is not a way around this either. Consent given for
research participation is not consent to be a test fixture in an unrelated
project, and openSNP's 2025 shutdown was driven in part by exactly this kind of
downstream reuse.

## Writing a fixture

Fixtures are hand-written to exercise a specific parsing behaviour. Keep them
small — a dozen lines is usually plenty — and make each one obvious about what
it is testing.

Use these reserved rsIDs for invented markers so a fixture can never collide
with a real variant in the catalog:

| Range | Purpose |
|---|---|
| `rs0`–`rs99` | Invented markers. dbSNP does not issue rsIDs in this range. |
| `i0`–`i99` | Invented 23andMe-style internal identifiers. |

Real rsIDs are fine when a test genuinely needs one — testing that APOE
derivation works requires `rs429358` — because the genotypes attached to them
are invented. What makes a file safe is that no real person's calls are in it.

Naming: `<provider>_<what_it_tests>.txt`, for example `23andme_no_calls.txt`,
`ancestry_numeric_chromosomes.txt`.

## Checking before you commit

The repo `.gitignore` blocks `*.txt` and `*.csv` by default and un-ignores this
directory specifically, so a stray export elsewhere in the tree will not be
staged by accident. That protects against the common mistake, not a deliberate
`git add -f`.

Before committing a fixture:

```sh
# Real exports are hundreds of thousands of lines. A fixture is dozens.
wc -l tests/fixtures/*

# Whatever you are adding should be something you can read in full.
cat tests/fixtures/your_new_fixture.txt
```

If a file is large enough that you cannot eyeball it, that is the signal: it is
not a fixture.
