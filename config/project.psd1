@{
    ProjectName    = "CHANGE_ME"
    AnalysisProfile = "Primary_padj0.05_LFC1.0"
    CountMatrix    = "INPUT/matrix_gene.count.xls"
    Metadata       = "INPUT/sample_metadata.csv"
    WorkRoot       = "work"
    ResultsRoot    = "results"
    SharedCache    = "resources/cache"
    DefaultPadjThreshold = 0.05
    DefaultLfcThreshold  = 1.0

    ControlGroup   = "Control"
    TreatGroup     = "Treatment"
    ExcludeSamples = ""
    BatchColumn    = ""
    MinCount       = 10
    MinSamples     = 3
    HeatmapTopPerDirection = 20

    PrepareCache   = $false
    EnablePathview = $false
}
