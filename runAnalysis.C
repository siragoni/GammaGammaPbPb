// to run:
//   alien-token-init
//   source /tmp/gclient_env_501
//   aliroot runAnalysis.C\(opt\)


// include the header of your analysis task here! for classes already compiled by aliBuild,
// precompiled header files (with extension pcm) are available, so that you do not need to
// specify includes for those. for your own task however, you (probably) have not generated a
// pcm file, so we need to include it explicitly
#include "AliAnalysisTaskNanoJPsi2016Fwd.h"

void runAnalysis(Int_t opt)
{


  /**
   * RUN ON LOCAL:
   */
  // Bool_t local    = kTRUE;
  // Bool_t gridTest = kFALSE;

  /**
   * RUN ON GRIDTEST:
   */
  // Bool_t local    = kFALSE;
  // Bool_t gridTest = kTRUE;


  /**
   * FULL GRID MOD:
   */
  Bool_t local    = kFALSE;
  Bool_t gridTest = kFALSE;





    // since we will compile a class, tell root where to look for headers
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    // compile the class and load the add task macro
    // here we have to differentiate between using the just-in-time compiler
    // from root6, or the interpreter of root5
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("AliAnalysisTaskNanoJPsi2016Fwd.cxx++g");
    AliAnalysisTaskNanoJPsi2016Fwd *task = reinterpret_cast<AliAnalysisTaskNanoJPsi2016Fwd*>(gInterpreter->ExecuteMacro("AddNanoJPsi2016Fwd.C"));
#else
    gROOT->LoadMacro("AliAnalysisTaskNanoJPsi2016Fwd.cxx++g");
    gROOT->LoadMacro("AddNanoJPsi2016Fwd.C");
    AliAnalysisTaskUPCforwardMC *task = AddTaskUPCforwardMC();
#endif


    if(!mgr->InitAnalysis()) return;
    // mgr->SetDebugLevel(2);
    // mgr->PrintStatus();
    // mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        // add a few files to the chain (change this so that your local files are added)

        // FILE *f = fopen("listMCsample.txt","r");
        FILE *f = fopen("listAxE.txt","r");
        // FILE *f = fopen("list.txt","r");
        char fileadd[300];
        Int_t flaggingValue = 0;
        while(fscanf(f,"%s",fileadd)==1){
            // chain->AddFile(fileadd);
            flaggingValue = chain->Add(fileadd);
            if(flaggingValue == 0) std::cout << fileadd << std::endl;
            flaggingValue = 0;
        }


        // chain->Add("AliAOD/*");

        // chain->Add("LHC18q/0001/*");
        // chain->Add("LHC18q/0002/*");
        // chain->Add("LHC18q/0003/*");
        // chain->Add("LHC18q/0004/*");
        // chain->Add("LHC18q/0005/*");
        // chain->Add("LHC18q/0006/*");
        // chain->Add("LHC18q/0007/*");
        // chain->Add("LHC18q/0008/*");
        // chain->Add("LHC18q/0009/*");

        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskNanoJPsi2016Fwd.cxx AliAnalysisTaskNanoJPsi2016Fwd.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskNanoJPsi2016Fwd.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20181028_ROOT6-1");
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");
        // select the input data for 2018 q
        // alienHandler->SetCheckCopy(kFALSE);




        // select the input data for 2018 q
      	if (opt == 0) {
      	  alienHandler->SetGridDataDir("/alice/data/2018/LHC18q");
      	  alienHandler->SetDataPattern("*muon_calo_pass3/AOD225/PWGUD/UD_PbPb_AOD/501_20190723-1440/*AliAOD.UPCNano.root");
      	  // MC has no prefix, data has prefix 000
      	  alienHandler->SetRunPrefix("000");
      	  // runnumber
      	  alienHandler->AddRunNumber(296510);
          // for( Int_t iRunLHC18q = 0; iRunLHC18q < 128; iRunLHC18q++){
          //   alienHandler->AddRunNumber( listOfGoodRunNumbersLHC18q[iRunLHC18q] );
          // }

      	} else if (opt == 1) {
      	  alienHandler->SetGridDataDir("/alice/data/2018/LHC18r");
      	  alienHandler->SetDataPattern("*muon_calo_pass3/AOD225/PWGUD/UD_PbPb_AOD/502_20190723-1441/*AliAOD.UPCNano.root");
      	  // MC has no prefix, data has prefix 000
      	  alienHandler->SetRunPrefix("000");
      	  // runnumber
      	  alienHandler->AddRunNumber(297481);
          // for( Int_t iRunLHC18r = 0; iRunLHC18r <  97; iRunLHC18r++){
          //   alienHandler->AddRunNumber( listOfGoodRunNumbersLHC18r[iRunLHC18r] );
          // }

      	} else if (opt == 2) {
      	  alienHandler->SetGridDataDir("/alice/data/2015/LHC15o");
      	  alienHandler->SetDataPattern("*muon_calo_pass1/AOD229/PWGUD/UD_PbPb_AOD/522_20190828-1630/*AliAOD.UPCNano.root");
      	  // MC has no prefix, data has prefix 000
      	  alienHandler->SetRunPrefix("000");
      	  // runnumber
      	  alienHandler->AddRunNumber(246087);
          // for( Int_t iRunLHC18r = 0; iRunLHC18r <  97; iRunLHC18r++){
          //   alienHandler->AddRunNumber( listOfGoodRunNumbersLHC18r[iRunLHC18r] );
          // }

      	} else {
      	  cout << " not a valid option ... bye!" << endl;
      	}


        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(50);
        alienHandler->SetExecutable("MC.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(10000);
        alienHandler->SetJDLName("MC.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        // merging: run with "kTRUE" and "full" for normal run
        // to merge on grid run jobs in SetRunMode("terminate")
        // to collect final results set SetMergeViaJDL(kFALSE)
        // alienHandler->SetMergeViaJDL(kTRUE);

        /* - The setting to kFALSE is to download the output files
           -
         */
        alienHandler->SetMergeViaJDL(kFALSE);
        alienHandler->SetMaxMergeStages(3);


        TString LHC18q("LHC18q");
        TString LHC18r("LHC18r");
        TString LHC15o("LHC15o");
        // define the output folders
        // alienHandler->SetGridWorkingDir("myWorkingDir");
        // alienHandler->SetGridWorkingDir("CheckAD_null");
        // alienHandler->SetGridWorkingDir("AD_HWandSW");
        // alienHandler->SetGridWorkingDir("PbPb15o_polPt2_test");
        if (opt == 0) alienHandler->SetGridWorkingDir("PbPb18q_pol_test7");
        if (opt == 1) alienHandler->SetGridWorkingDir("PbPb18r_pol_test7");
        if (opt == 2) alienHandler->SetGridWorkingDir("PbPb15o_pol_test7");

        // alienHandler->SetGridOutputDir("myOutputDir");
        if (opt == 0) alienHandler->SetGridOutputDir(LHC18q.Data());
        if (opt == 1) alienHandler->SetGridOutputDir(LHC18r.Data());
        if (opt == 2) alienHandler->SetGridOutputDir(LHC15o.Data());



        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis

            /* - The option FULL is to send the full analysis.
               -
             */
            // alienHandler->SetRunMode("full");

            /* - This option TERMINATE is used for the merging of the files.
               -
             */
            alienHandler->SetRunMode("terminate");
            mgr->StartAnalysis("grid");
        }
    }
}

/* -
 * - Welcome my dear ALICE user! To use ALICE software from CVMFS:
 * - List all packages         --> alienv q
 * - List AliPhysics packages  --> alienv q | grep -i aliphysics
 * - Enable a specific package --> alienv enter VO_ALICE@AliPhysics::vAN-20190114_ROOT6-1
 * - Enable multiple packages  --> alienv enter VO_ALICE@AliPhysics::vAN-20190114_ROOT6-1,VO_ALICE@fastjet::v3.2.1_1.024-alice3-7
 */
