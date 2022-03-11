#include "DataEvent.h"

DataEvent::DataEvent(string fileName, string treeName){
    SetUpTree(fileName, treeName);
    numofRootFiles=1;
}

DataEvent::DataEvent(string filesfileName, string treeName, int numOfFiles){
    numofRootFiles=numOfFiles;
    SetUpChain(filesfileName, treeName);
}

void DataEvent::SetUpTree(string fileName, string treeName){
    loc_Tree = new TChain(treeName.c_str());
    loc_Tree->Add(fileName.c_str());
    //cout<<"Setting up tree "<<treeName.c_str()<<" in root file "<<fileName.c_str()<<endl;
    SetUpBranches();
}

void DataEvent::SetUpChain(string filesfileName, string treeName){
    loc_Tree = new TChain(treeName.c_str());
    string loc_Name_File;
    ifstream myfile;
    myfile.open(filesfileName.c_str());
    loc_Name_File.clear();
    cout<<"adding files to chain:"<<endl;
    numofRootFiles=0;
    while (myfile >> loc_Name_File){
        numofRootFiles++;
        std::cout<<loc_Name_File.c_str()<<endl;
        loc_Tree->Add(loc_Name_File.c_str());
        loc_Name_File.clear();
    }
    cout<<"Added "<<numofRootFiles<<" files"<<endl;
    SetUpBranches();
    myfile.close();
}

void DataEvent::SetUpBranches(){
    loc_EVNT_track=0;
    loc_eloss_track=0;
    loc_kpi_mm=0;
    loc_k_mm=0;
    loc_ppi_mm=0;
    loc_pipi_mm=0;

    loc_ECPB_ECtot=0;
    loc_ECPB_ECinn=0;
    loc_ECPB_ECout=0;
    loc_ECPB_X=0;
    loc_ECPB_Y=0;
    loc_ECPB_Z=0;
    loc_ECPB_time=0;
    loc_ECPB_path=0;
    
    loc_EVNT_q=0;
    loc_EVNT_scsec=0;
    loc_EVNT_scpad=0;
    loc_EVNT_schit=0;
    loc_EVNT_stsec=0;
    loc_EVNT_sthit=0;
    loc_tagr_eid=0;
    loc_tagr_tid=0;
    loc_tagr_stat=0;
    loc_STPB_sthid=0;
    loc_SCPB_ScPdHt=0;
    loc_index_k=0;
    loc_index_pi=0;
    loc_EVNT_bem=0;
    loc_EVNT_sc_t=0;
    loc_EVNT_sc_d=0;
    loc_EVNT_st_t=0;
    loc_EVNT_st_d=0;
    loc_EVNT_sc_e=0;
    loc_delt_t_k=0;
    loc_delt_t_pi=0;
    loc_tagr_epho=0;
    loc_tagr_tpho=0;
    loc_coh_plan=0;
    loc_coh_radi=0;
    loc_coh_plan_db=0;
    loc_coh_edge=0;
    loc_beam_en=0;
    loc_coh_edge_nom=0;
    loc_numph_k=0;
    loc_numph_pi=0;

    loc_EVNT_vertex=0;
    loc_MVRT_vertex=0;
    eventno=0;
//    cout<<"Setting up branches"<<endl;
    loc_Tree->SetBranchAddress("EVNT_track",&loc_EVNT_track); //vector LorentzVector
    loc_Tree->SetBranchAddress("eloss_track",&loc_eloss_track); //vector Eloss corrected LorentzVector
    loc_Tree->SetBranchAddress("EVNT_q",&loc_EVNT_q);
    loc_Tree->SetBranchAddress("EVNT_scsec",&loc_EVNT_scsec);
    loc_Tree->SetBranchAddress("EVNT_scpad",&loc_EVNT_scpad);
    loc_Tree->SetBranchAddress("EVNT_schit",&loc_EVNT_schit);
    loc_Tree->SetBranchAddress("EVNT_stsec",&loc_EVNT_stsec);
    loc_Tree->SetBranchAddress("EVNT_sthit",&loc_EVNT_sthit);
    loc_Tree->SetBranchAddress("tagr_tid",&loc_tagr_tid);
    loc_Tree->SetBranchAddress("tagr_eid",&loc_tagr_eid);
    loc_Tree->SetBranchAddress("tagr_stat",&loc_tagr_stat);
    loc_Tree->SetBranchAddress("SCPB_ScPdHt",&loc_SCPB_ScPdHt);
    loc_Tree->SetBranchAddress("STPB_sthid",&loc_STPB_sthid);
    loc_Tree->SetBranchAddress("EVNT_bem",&loc_EVNT_bem);
    loc_Tree->SetBranchAddress("EVNT_sc_t",&loc_EVNT_sc_t);
    loc_Tree->SetBranchAddress("EVNT_sc_d",&loc_EVNT_sc_d);
    loc_Tree->SetBranchAddress("EVNT_st_t",&loc_EVNT_st_t);
    loc_Tree->SetBranchAddress("EVNT_st_d",&loc_EVNT_st_d);
    loc_Tree->SetBranchAddress("EVNT_sc_e",&loc_EVNT_sc_e);
    loc_Tree->SetBranchAddress("tagr_epho",&loc_tagr_epho);
    loc_Tree->SetBranchAddress("tagr_tpho",&loc_tagr_tpho);
    loc_Tree->SetBranchAddress("numofpart", &loc_numofpart);
    loc_Tree->SetBranchAddress("num_pos", &loc_num_pos);
    loc_Tree->SetBranchAddress("num_neg", &loc_num_neg);
    loc_Tree->SetBranchAddress("num_neu", &loc_num_neu);
    loc_Tree->SetBranchAddress("num_deuterons", &loc_num_deuterons);
    loc_Tree->SetBranchAddress("num_protons", &loc_num_protons);
    loc_Tree->SetBranchAddress("num_poskaons", &loc_num_poskaons);
    loc_Tree->SetBranchAddress("num_pospions", &loc_num_pospions);
    loc_Tree->SetBranchAddress("num_negkaons", &loc_num_negkaons);
    loc_Tree->SetBranchAddress("num_negpions", &loc_num_negpions);
    loc_Tree->SetBranchAddress("head_eventnum", &loc_head_eventnum); //Event number
    loc_Tree->SetBranchAddress("head_runnum", &loc_head_runnum); //Run Number
    loc_Tree->SetBranchAddress("trip_flag", &loc_trip_flag);
    loc_Tree->SetBranchAddress("EVNT_vertex",&loc_EVNT_vertex);
    loc_Tree->SetBranchAddress("MVRT_vertex",&loc_MVRT_vertex);
    loc_Tree->SetBranchAddress("coh_edge",&loc_coh_edge);
    loc_Tree->SetBranchAddress("beam_en",&loc_beam_en);
    loc_Tree->SetBranchAddress("coh_edge_nom",&loc_coh_edge_nom);
    loc_Tree->SetBranchAddress("coh_plan_db",&loc_coh_plan_db);
    loc_Tree->SetBranchAddress("coh_radi",&loc_coh_radi);
    loc_Tree->SetBranchAddress("coh_plan",&loc_coh_plan);
    loc_Tree->SetBranchAddress("delt_t_k",&loc_delt_t_k);
    loc_Tree->SetBranchAddress("delt_t_pi",&loc_delt_t_pi);
    loc_Tree->SetBranchAddress("numph_k",&loc_numph_k);
    loc_Tree->SetBranchAddress("numph_pi",&loc_numph_pi);
    loc_Tree->SetBranchAddress("index_k",&loc_index_k);
    loc_Tree->SetBranchAddress("index_pi",&loc_index_pi);
    loc_Tree->SetBranchAddress("kpi_mm",&loc_kpi_mm);
    loc_Tree->SetBranchAddress("k_mm",&loc_k_mm);
    loc_Tree->SetBranchAddress("ppi_mm",&loc_ppi_mm);
    loc_Tree->SetBranchAddress("pipi_mm",&loc_pipi_mm);
    loc_Tree->SetBranchAddress("ECPB_ECtot",&loc_ECPB_ECtot);
    loc_Tree->SetBranchAddress("ECPB_ECinn",&loc_ECPB_ECinn);
    loc_Tree->SetBranchAddress("ECPB_ECout",&loc_ECPB_ECout);
    loc_Tree->SetBranchAddress("ECPB_X",&loc_ECPB_X);
    loc_Tree->SetBranchAddress("ECPB_Y",&loc_ECPB_Y);
    loc_Tree->SetBranchAddress("ECPB_Z",&loc_ECPB_Z);
    loc_Tree->SetBranchAddress("ECPB_time",&loc_ECPB_time);
    loc_Tree->SetBranchAddress("ECPB_path",&loc_ECPB_path);
    
//    cout<<"Branches setup"<<endl;
    loc_Tree->GetEntries();
//    cout<<"File has "<<loc_Tree->GetEntries()<<" events"<<endl;
}

int DataEvent::getNextEntry(){
    if (eventno<loc_Tree->GetEntries()){
        loc_Tree->GetEvent(eventno);
        eventno++;
        return eventno;
    }else return -1;
}
int DataEvent::getEntries(){
    return loc_Tree->GetEntries();
}


float DataEvent::getECPB_pathcorr(int ip){
    float ECtotal=getECPB_ECtot(ip);
    float ECouter=getECPB_ECout(ip);
    float ECinner=getECPB_ECinn(ip);
    float newpath=0;
//    if( strcmp(opt,”g13a”) == 0 ){
//            //EC average corrections obtained from gd->pi+pi-p(n) analysis
//        if((ECtotal>0.0)&&(ECinner>0.0)&&(ECouter>0.0))
//            ecpath_neutral=ecpath_neutral+15.1986;//g13a
//        if((ECtotal>0.0)&&(ECinner>0.0)&&(ECouter==0.0))
//            ecpath_neutral=ecpath_neutral+6.09480;//g13a
//        if((ECtotal>0.0)&&(ECinner==0.0)&&(ECouter>0.0))
//            ecpath_neutral=ecpath_neutral-8.42826;//g13a
//    }
//    else if( strcmp(opt,”g13b”) == 0 ){
//            //EC average corrections obtained from gd->pi+pi-p(n) analysis
    if((ECtotal>0.0)&&(ECinner>0.0)&&(ECouter>0.0))
        newpath = getECPB_path(ip) + 17.86;//g13b
    if((ECtotal>0.0)&&(ECinner>0.0)&&(ECouter==0.0))
        newpath = getECPB_path(ip) + 9.78; //g13b
    if((ECtotal>0.0)&&(ECinner==0.0)&&(ECouter>0.0))
        newpath = getECPB_path(ip) - 8.76; //g13b
//    }
//    else if( strcmp(opt,”MC”) == 0 ){
//            //EC average corrections obtained from gd->pi+pi-p(n) analysis
//        if((ECtotal>0.0)&&(ECinner>0.0)&&(ECouter>0.0))
//            ecpath_neutral = ecpath_neutral + 55.5114;//MC
//        if((ECtotal>0.0)&&(ECinner>0.0)&&(ECouter==0.0))
//            ecpath_neutral = ecpath_neutral + 51.9899;//MC
//        if((ECtotal>0.0)&&(ECinner==0.0)&&(ECouter>0.0))
//            ecpath_neutral = ecpath_neutral + 26.6853;//MC
//    }
//    else{
//        cout << “Wrong option in NeutronCorrections.C: check <opt> argument” << endl;
//        exit(1);
//    }
    return newpath;
}

float DataEvent::getECPB_betacorr(int ip){
    double newbeta=0;
    double startime=0;
    double speedOflight=29.9792458;
    startime=getEVNT_sc_t(1)-getEVNT_sc_d(1)/(getEVNT_bem(1)*speedOflight);
    newbeta=getECPB_pathcorr(ip)/(getECPB_time(ip)-startime)/speedOflight;
    return newbeta;
}



