if VERSION >= v"1.8" # problems with previous Julia versions...
    @setup_workload begin
        # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
        # precompile file and potentially make loading faster.
        gI = [1; 3; 5;; 2; 4; 5;;; 1; 4; 5;; 2; 5; 5];
        bI = [2; 4; 2;; 3; 7; 5;;; 2; 3; 2;; 3; 6; 5];
        gO = [10; 30; 50;; 20; 40; 50;; 15; 8; 12;;; 12; 40; 50;; 22; 50; 50;; 16; 55; 55];
        bO = [2; 4; 2;; 3; 7; 5;;; 2; 3; 2;; 3; 6; 5];
        X = [4 7 8 4 2 10; 3 3 1 2 4 1]
        Y = [1 1 1 1 1 1]
        @compile_workload begin
            @info "Beginning BDisposal PrecompileTool workflow...."
            prodIndexFB(gI,gO,bO,bI;remarcable_obs_dmu=1, remarcable_obs_period=2);
            prodIndex(gI,gO,bO,bI,retToScale="variable",convexAssumption=false)
            prodIndex(gI,gO,bO,bI,retToScale="variable",convexAssumption=true,prodStructure="additive")
            prodIndex(gI,gO,bO,bI,retToScale="variable",convexAssumption=true,prodStructure="multiplicative")
            efficiencyScores(gI,gO,bO,bI,retToScale="variable", dirGI=0,dirBI=0,dirGO=1,dirBO=0, prodStructure="multiplicative")
            dmuEfficiency(X[:,1],Y[:,1],X',Y')
            dmuEfficiencyDual(X[:,1],Y[:,1],X',Y')
            @info "...done BDisposal PrecompileTool workflow."
        end
    end
end