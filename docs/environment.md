# Environment assessment

Assessment date: 2026-11-27. This inventory records availability, not scientific-environment qualification.

## Observed baseline

| Component | Observation | Qualification still needed |
| --- | --- | --- |
| Host | macOS / Darwin, arm64 | Compute resources and execution constraints for real cases |
| Git | 2.55.0; repository already initialized on `main`, initially without commits or remotes | None for initialization |
| GitHub CLI | 2.102.0; authenticated account `K1000-B`; HTTPS Git protocol | Network access may require a narrowly scoped execution permission |
| Docker CLI | 29.8.0, build `88096ef` | Daemon access, intended container/image, SU2 version/build, MPI and launch conventions |
| Gmsh | Executable found | Version and required native-language capabilities |
| ParaView | Executable found | Version and intended reproducible visualization workflow |
| Python | `python3` executable found | User-designated environment, interpreter version and approved dependencies |
| SU2 | `SU2_CFD` not found on the inspected host PATH | Intended Docker-based solver has not been supplied or inspected |

No Python interpreter was launched for initialization. No environment was created, package installed, application reconfigured, Docker daemon queried, mesh generated or CFD simulation launched. Executable presence alone does not confirm runtime readiness.

The first sandboxed `gh auth status` reported an invalid token; repeating the check with network-capable execution succeeded, and `gh api user` confirmed the account. Do not request a new login merely from that restricted check if the network-capable check succeeds.

## Available collaboration mechanisms

The current Codex session exposes filesystem/terminal tools, official-web lookup, connected-app tools and desktop/browser capabilities. Tool availability and permissions may differ in later sessions; inspect them rather than assuming persistence. The GitHub connector exposes repository/issue/PR operations but no repository-creation operation was found in the inspected tool inventory; the authenticated CLI handles initial repository creation and Git transport.

Durable project continuity uses repository files and Git. The root `AGENTS.md` directs future sessions to the current state and relevant documents. This follows the official [AGENTS.md discovery mechanism](https://learn.chatgpt.com/docs/agent-configuration/agents-md). No global settings, custom skills, hooks, scheduled tasks or project-specific Codex configuration are needed for the current scope.

## Information required before scientific execution

1. User-designated Python environment and permitted dependency-management procedure.
2. Intended SU2 Docker image/container and exact solver version; approved launch/resource constraints. Inspect only with the appropriate access and do not modify the environment.
3. Gmsh and ParaView versions when their features are first used.
4. Existing adaptation code and supporting documents, with provenance and publication authorization.
5. Candidate airfoil and experimental references, followed by justified physical and numerical choices.

Check SU2 options against the actual version's official template/source and Gmsh syntax against official documentation when implementation begins. No solver or library parameters are currently prescribed.
