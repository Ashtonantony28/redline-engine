# CLAUDE.md: conventions for Redline

1. **Every change ships with a Verify command** in its PR description (pytest invocation or module run).
2. **Analytic before numeric.** Each physics model has an analytic or published reference test in `tests/`
   before it is used downstream. Tolerances are stated, never loosened silently.
3. **`spec/engine.yaml` is the single source of truth**; no magic numbers in model code.
4. Units and conventions follow ADR 0001. Hard-to-reverse decisions get an ADR in `docs/adr/`.
5. All results are **modelled**: there is no physical engine. Say so wherever numbers are reported.
6. Never commit to `main` directly once the GitHub repo exists; one logical change per PR; red CI never
   merges; commit messages explain why.
