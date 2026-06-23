@{
    ProjectName    = "CHANGE_ME"
    AnalysisProfile = "Exploratory_padj0.05_LFC0.5"
    CountMatrix    = "INPUT/matrix_gene.count.xls"
    Metadata       = "INPUT/sample_metadata.csv"
    WorkRoot       = "work"
    ResultsRoot    = "results"
    SharedCache    = "resources/cache"
    DefaultPadjThreshold = 0.05
    DefaultLfcThreshold  = 0.5

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
