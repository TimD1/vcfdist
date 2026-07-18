### Creating Issues

### Creating Feature Branches
 - All work should be completed on a feature branch in a git worktree, and squash-merged into dev by the user. DO NOT merge work yourself.
 - Feature branches should be named "<issue>_td_<brief-description>", where <brief-description> will usually be prefixed with "D<deliverable_id>-"

### Documentation Conventions
 - .h files: /** @brief ... */ only — no @param, @return, @throws
 - .cpp files: full /** @brief ... @param[in] ... @return ... @throws ... */
 - @param[in] for all input parameters; @param[out] for output-only; @param[in,out] for modified
 - @return (not @returns)
 - @throws ERROR <description> whenever the function calls ERROR()
 - @throws WARNING <description> whenever the function calls WARN()
 - Separator lines: exactly 100 characters wide (e.g. /* + 96 * + /, or /* Label ***...*/ padded to 100)
