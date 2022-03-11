enum {
    E_ID    = 0,
    ENERGY  = 1,
    ENH     = 2,
    ENHERR  = 3,
    ENHFIT  = 4,
    PFIT    = 5,
    PCOR    = 6,
    PCORERR = 7,
    PSMOOTH = 8,
    EDGE    = 9
};

enum {
    PARA1=0,
    PERP1=1,
    PARA2=2,
    PERP2=3,
    PARA3=4,
    PERP3=5,
    PARA4=6,
    PERP4=7,
    PARA5=8,
    PERP5=9,
    PARA6=10,
    PERP6=11,
    PARA7=12,
    PERP7=13,
    PARA8=14,
    PERP8=15,
    PARA9=16,
    PERP9=17,
    PARA10=18,
    PERP10=19,
    PARA11=20,
    PERP11=21,
    PARA12=22,
    PERP12=23};

Double_t polTable[24][500][385][10];  //where its [plane][edge][E_id][field]
Int_t polTableN[24]={0,0,0,0,0,0};            //No of entries for para and perp
Char_t polFirstLines[24][500][250];   //to keep the 1st lines if needed (info from files)

Int_t edgeEventLow[5000];            //hold the current table of edge positions for event ranges
Int_t edgeEventHigh[5000];
Double_t edgeEventEdge[5000];
Int_t edgeEventPlane[5000];
Int_t edgeEventN;
Int_t edgeIndex=0;
Int_t lastEdgeEvent=0;
Double_t lastCohEdge=0.0;
Int_t lastCohPlane=-1;



int LoadPolTable(int plane, Char_t *PolTableList){
    FILE *fplist,*fpfile;              //file pointers for filelist and each file therein
    Char_t lline[250];                 //for reading in lines from filelist
    Char_t fline[250];                 //for reading in lines from file
    Char_t filename[250];              //file
    Int_t  fcount=0;                   //counter for no of files read in
    Int_t  chancount=0;                //counter for no of channels read in
    Int_t  eid=0;
    Double_t edge=0.0;

    if((fplist=fopen(PolTableList,"r"))==NULL){ //open filelist
        cerr << "Error Couldn't open file: " << PolTableList << endl;
        return -1;
    }

    fcount=0;
    //for each file in the list
    while(fgets(lline,240,fplist) != NULL){
        if((lline[0] == '*')||(lline[0] == '#')) continue; //skip comments
        sscanf(lline,"%s",filename);                       //read in filename

        cout << "opening " << filename << "   " << endl;
        if((fpfile=fopen(filename,"r"))==NULL){              //open file
            cerr << "Error Couldn't open file: " << filename << endl;
            return -1;
        }

        fgets(polFirstLines[plane][polTableN[plane]],240,fpfile ); //save the 1st line

        //scan the bit after the last "_" in the filename to get the edge energy
        sscanf(strrchr(filename,'_')+1,"%lg",&polTable[plane][fcount][0][EDGE]);

        chancount=0;                                             //starting array at 1 for consistency with E_ID
        while((fgets(fline,240,fpfile)) != NULL){
            if((fline[0] == '*')||(fline[0] == '#')) continue;     //skip comments
            sscanf(fline,"%d",&eid);                               //first get the E_ID
            sscanf(fline,"%*d%lg%lg%lg%lg%lg%lg%lg%lg",
                   &polTable[plane][fcount][eid][ENERGY],
                   &polTable[plane][fcount][eid][ENH],
                   &polTable[plane][fcount][eid][ENHERR],
                   &polTable[plane][fcount][eid][ENHFIT],
                   &polTable[plane][fcount][eid][PFIT],
                   &polTable[plane][fcount][eid][PCOR],
                   &polTable[plane][fcount][eid][PCORERR],
                   &polTable[plane][fcount][eid][PSMOOTH]);
            chancount++;
        }
        fclose(fpfile); //close the file
        if(chancount!=384){
            cerr << "Should be 384 lines in " << filename << " - only got " << chancount << endl;
            return -1;
        }
        polTableN[plane]++;

        fcount++;
    }

    fclose(fplist);

    return(0);
}


Double_t GetPol(Int_t plane, Double_t edge, Int_t eid, Int_t poltype = PSMOOTH, Double_t lowThresh=0.2, Double_t highThresh=0.3){
    //get polarization based on eid and edge position

    Int_t eIndex=0;
    Double_t pol=-1.0;

    //Check edge in in range of tables
    if((edge<polTable[plane][1][0][EDGE])||(edge>polTable[plane][polTableN[plane]-1][0][EDGE])) return -1.0;

    //cout << "In range" << endl;

    //find index
    for(eIndex=0;eIndex<polTableN[plane];eIndex++){
        if(polTable[plane][eIndex][0][EDGE]>=edge) break;
    }
    //cout << "Index = " << eIndex << endl;

    pol=polTable[plane][eIndex][eid][poltype];
    if((polTable[plane][eIndex][0][ENERGY]<edge)&&(pol<lowThresh)) pol = -1.0;
    if((polTable[plane][eIndex][0][ENERGY]>edge)&&(pol<highThresh)) pol = -1.0;

    return pol;
}


Double_t GetPol(Int_t plane, Double_t edge, Double_t eg, Int_t poltype = PSMOOTH, Double_t lowThresh=0.2, Double_t highThresh=0.3){
    //get polarization based on ephoton energy and edge position

    Int_t eIndex=0;
    Double_t pol=-1.0;
    Int_t eid=0;

    //Check edge in in range of tables
    if((edge<polTable[plane][1][0][EDGE])||(edge>polTable[plane][polTableN[plane]-1][0][EDGE])) return -1.0;

    //find index
    for(eIndex=0;eIndex<polTableN[plane];eIndex++){
        if(polTable[plane][eIndex][0][EDGE]>=edge) break;
    }
    //cout << "Index = " << eIndex << endl;

    //find eid
    for(eid=1;eid<=384;eid++){
        if(polTable[plane][eIndex][eid][ENERGY]<=eg) break;
    }
     if (eg >= (polTable[plane][eIndex][eid][ENERGY]+polTable[plane][eIndex][eid-1][ENERGY])/2.0)
        eid=eid-1;


    //cout << "eid = " << eid <<endl;

    pol=polTable[plane][eIndex][eid][poltype];
    if((polTable[plane][eIndex][0][ENERGY]<edge)&&(pol<lowThresh)) pol = -1.0;
    if((polTable[plane][eIndex][0][ENERGY]>edge)&&(pol<highThresh)) pol = -1.0;

    return pol;
}
