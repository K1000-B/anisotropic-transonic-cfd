# Project status and decisions

Last updated: 2026-11-27.

## Current state

Initialization work only. Permanent instructions, scientific scope, workflow and an environment assessment are documented. No scientific implementation, case, mesh, imported dataset or computed result is present. No scientific claim has been tested or physically validated.

GitHub target: `K1000-B/anisotropic-transonic-cfd`, intended public URL: https://github.com/K1000-B/anisotropic-transonic-cfd. Creation and initial synchronization are pending at this point; update this status after verification. Git history is the authoritative record of actual commits.

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
- Content review, link/whitespace checks, initial commit, remote creation and publication verification remain to be completed before marking initialization finished.

## Next actions

1. Complete and verify initial Git/GitHub publication, then record the outcome here.
2. Receive the existing algorithm and supporting materials outside the public repository, clarify provenance and publication rights, and perform a read-only technical assessment before a scoped import.
3. Receive the designated Python environment and SU2 Docker details. Do not create or modify them in anticipation.
4. Use the actual algorithm and primary experimental sources to define the next scientific task and arbitrate important physical choices with the user.

## Open issues and scientific cautions

- Publication rights and licensing of external work are unknown.
- Scientific software versions, Python dependencies, computational resources and exact SU2 launch workflow are not yet qualified.
- Airfoil, reference experiment, physical conditions and numerical settings remain undecided.
- Exterior Euler-mesh reuse for RANS must be checked for viscous adequacy; wall discretization and shock displacement may require targeted changes.
- Future convergence and shock-extraction methods must account for discontinuities and nonuniform anisotropic mesh sequences; apparent convergence alone is not validation.
