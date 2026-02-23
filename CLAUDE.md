

### Documentation Conventions
 - .h files: /** @brief ... */ only — no @param, @return, @throws
 - .cpp files: full /** @brief ... @param[in] ... @return ... @throws ... */
 - @param[in] for all input parameters; @param[out] for output-only; @param[in,out] for modified
 - @return (not @returns)
 - @throws ERROR <description> whenever the function calls ERROR()
 - @throws WARNING <description> whenever the function calls WARN()
 - Separator lines: exactly 100 characters wide (e.g. /* + 96 * + /, or /* Label ***...*/ padded to 100)
