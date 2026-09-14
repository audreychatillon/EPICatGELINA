void savetcutg_A03_v82(){

    TFile * fsave = new TFile("fsave_tcutg_discri_A03.root","recreate");
    for(int c = 0 ; c < 4; c++){
        gInterpreter->ExecuteMacro(Form("A03_cut%02d.C", c+1));
        TCutG * cut = (TCutG*)gROOT->FindObject(Form("A03_cut%02d",c+1));
        fsave->cd();
        cut->Write();
    }
    fsave->ls();
    fsave->Close();

}
