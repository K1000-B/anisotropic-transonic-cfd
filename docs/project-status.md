# Project status and decisions

Last updated: 2026-11-27.

## Current state

Initialization work only. Permanent instructions, scientific scope, workflow and an environment assessment are documented. No scientific implementation, case, mesh, imported dataset or computed result is present. No scientific claim has been tested or physically validated.

Public repository: [K1000-B/anisotropic-transonic-cfd](https://github.com/K1000-B/anisotropic-transonic-cfd). Created and initially synchronized on 2026-11-27. `origin` uses HTTPS and `main` tracks `origin/main`. Description and seven scientific topics are configured. Git history is the authoritative record of actual commits.

## Organization decisions

| Decision | Rationale |
| --- | --- |
| Root `AGENTS.md` plus focused documentation | Durable session entry and clear ownership of rules, science, environment and progress without relying on chat memory |
| Documentation-only initial tree | Existing scientific work has not arrived; avoid speculative source/case directories and empty scaffolding |
| `main`, atomic Conventional Commits, ordinary pushes | Simple traceable workflow suited to the current project size; branches and review can be added when useful |
| Repository-local GitHub noreply identity | Preserve commit attribution without publishing a personal email or changing global configuration |
| Ignore generated root output locations, not entire scientific formats | Keep large reproducible results out of Git while allowing reviewed reference meshes and datasets |
| No custom skill, automation, CI or dependency specification | No demonstrated repeated workflow or executable implementation yet; avoid unnecessary infrastructure |
| No license applied | Ownership and redistribution rights of pre-existing work remain unresolved; explicit owner approval is required |

A proposed license for original project code is [MIT](https://choosealicense.com/licenses/mit/), for its short permissive terms. This is a proposal only. Before applying it, confirm ownership/coauthor permissions and whether third-party code imposes other terms. Documentation and datasets may require separate terms and attribution. This initialization does not approve publication of the external algorithm or teaching material.

## Initial checks

- Inspected the existing repository, branch, remotes and pre-existing untracked operating-system artifact; no prior scientific files or commits were present.
- Checked installed executable availability and Git/GitHub CLI/Docker CLI versions without modifying scientific environments.
- Verified GitHub CLI authentication with network-capable execution and confirmed the account through the GitHub API.
- Consulted official Codex instructions documentation for repository-level continuity.
- Reviewed all seven initial files and the staged diff; checked local Markdown links, final newlines, trailing whitespace, and common local-path/credential patterns. No issue was found in these bounded checks; they are not a comprehensive automated security audit.
- Checked ignore behavior: local scratch material, generated outputs and Python environments are ignored; representative geometry, solver configuration, reference data and mesh inputs remain eligible for tracking.
- Created and pushed the initial commit `f6f651b` (`chore: initialize scientific project foundations`). Verified a clean working tree and identical local/remote `main` commit IDs after fetching.
- Confirmed public visibility and default branch through the GitHub API; an unauthenticated API request independently returned the repository metadata with `private: false`. Verified published root entries and configured topics.
- This status update records the publication outcome in a separate documentation commit. No scientific tests or physical validation were performed because no scientific implementation or case exists yet.

## Next actions

1. Receive the existing algorithm and supporting materials outside the public repository, clarify provenance and publication rights, and perform a read-only technical assessment before a scoped import.
2. Receive the designated Python environment and SU2 Docker details. Do not create or modify them in anticipation.
3. Resolve the proposed license with the user before applying any license grant.
4. Use the actual algorithm and primary experimental sources to define the next scientific task and arbitrate important physical choices with the user.

## Open issues and scientific cautions

- Publication rights and licensing of external work are unknown.
- Scientific software versions, Python dependencies, computational resources and exact SU2 launch workflow are not yet qualified.
- Airfoil, reference experiment, physical conditions and numerical settings remain undecided.
- Exterior Euler-mesh reuse for RANS must be checked for viscous adequacy; wall discretization and shock displacement may require targeted changes.
- Future convergence and shock-extraction methods must account for discontinuities and nonuniform anisotropic mesh sequences; apparent convergence alone is not validation.
