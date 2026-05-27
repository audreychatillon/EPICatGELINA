#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <Riostream.h>
#include <TROOT.h>

void run(int nDets)
{
    ofstream file_gamma(Form("default_Gamma_Peak.cal"));
    ofstream file_alpha(Form("default_Alpha1Dcut.cal"));
    ofstream file_tcutg(Form("default_TCutG_pts.cal"));

    for(int det = 1 ; det <= nDets ; det++){
        for(int anode = 1 ; anode <= 11 ; anode++){
            file_gamma << "EPIC_" << det << "_ANODE_" << anode << "_GAMMA_PEAK   0." << endl;  
            file_alpha << "EPIC_" << det << "_ANODE_" << anode << "_ALPHA 10000." << endl;        
            file_tcutg << "EPIC_" << det << "_ANODE_" << anode << "_TCUTG_DISCRI_X 25000. 50000. 50000. 25000." << endl; 
            file_tcutg << "EPIC_" << det << "_ANODE_" << anode << "_TCUTG_DISCRI_Y 1.0    1.0    0.5    0.5" << endl; 
            file_tcutg << endl ;
        }
        file_gamma << endl;  
        file_alpha << endl;        
    }

    file_gamma.close();
    file_alpha.close();
    file_tcutg.close();
}
