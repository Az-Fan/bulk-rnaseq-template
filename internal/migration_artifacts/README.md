# Historical migration artifacts

Project-specific source snapshots stored below this directory are read-only evidence for migration audits. They are not imported by `workflow/`, are not runtime dependencies, and are not a second environment or package-management path.

The FKBP1A-PAI snapshot deliberately excludes `.r-library*`, `work/`, `results/`, raw inputs and resources. Its official results are restored separately under the local project and verified by SHA256.
