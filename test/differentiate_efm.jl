@testset "Implicit differentiation of EFMs" begin 
    #: Set problem up 
    m = StandardModel("EFMModel")
    ms = ["A","A_ext","B","B_ext","C","D","E","E_ext","P","P_ext"]
    add_metabolites!(m,[Metabolite])


end