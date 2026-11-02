# BellOptimizer 1.2.122

## Consolidated Chinese-bell audit

Version 1.2.122 is a project-cleanup release. It does not introduce a new Chinese-bell physical model.

### Changes

- Consolidated the Chinese two-tone investigations from 1.2.102–1.2.121 into `docs/CHINESE_TWO_TONE_CONSOLIDATED_AUDIT.md`.
- Added `docs/CHINESE_BELL_LITERATURE.md` as the maintained bibliography.
- Removed experiment-specific result documents and exploratory sweep executables that are now superseded by the consolidated conclusions.
- Retained focused numerical/contract tests that protect the actual implementation rather than historical experiment logs.
- Updated project status and development plan to reference the consolidated audit.
- Production Western objective and structural FEM implementation are unchanged.

### Scientific status

The Chinese two-tone mechanism remains **PROVISIONAL / EXPERIMENTAL**. The current reduced FEM model has not yet demonstrated a robust physical `(m,n)a/(m,n)b` modal family together with reciprocal A/B strike selectivity.

No new physical thresholds or optimizer parameters were introduced.
