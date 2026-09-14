void savetcutg_A04_v4b(){

    TFile * fsave = new TFile("fsave_tcutg_discri_A04.root","recreate");
    for(int c = 0 ; c < 4; c++){
        gInterpreter->ExecuteMacro(Form("A04_cut%02d.C", c+1));
        TCutG * cut = (TCutG*)gROOT->FindObject(Form("A04_cut%02d",c+1));
        fsave->cd();
        cut->Write();
    }
    fsave->ls();
    fsave->Close();

}
