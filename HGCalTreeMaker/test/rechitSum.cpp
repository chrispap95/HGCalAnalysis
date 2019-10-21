void rechitSum() {
    TCanvas* c = new TCanvas("c","c",1);
    TFile* f = TFile::Open("out250.root");
    TTree* t = (TTree*)f->Get("hgcalTupleTree/tree");

    std::vector<float> *rechit = 0;
    std::vector<float> *rechiteta = 0;
    std::vector<float> *rechitphi = 0;
    std::vector<float> *rechitz = 0;
    std::vector<int> *rechitlayer = 0;
    std::vector<float> *etagen = 0;
    std::vector<float> *phigen = 0;
    std::vector<float> *ptgen = 0;

    t->SetBranchAddress("HGCRecHitEnergy",&rechit);
    t->SetBranchAddress("HGCRecHitEta",&rechiteta);
    t->SetBranchAddress("HGCRecHitPhi",&rechitphi);
    t->SetBranchAddress("HGCRecHitLayer",&rechitlayer);
    t->SetBranchAddress("HGCRecHitPosz",&rechitz);
    t->SetBranchAddress("GenParEta",&etagen);
    t->SetBranchAddress("GenParPhi",&phigen);
    t->SetBranchAddress("GenParPt",&ptgen);

    TH1F* h1 = new TH1F("h1","sum of HGCRecHitEnergy (GenParEta > 1.5 && #Delta R<0.2 && layer < 35);rechitsum [GeV]",50,50,250);

    for (long i = 1; i <= t->GetEntries(); ++i) {
        t->GetEntry(i-1);
        float rechitsum = 0;
        for (int j = 0; j<rechit->size();j++){
            if(etagen->at(0) > 1.5 && rechitlayer->at(j) < 29 && rechitz->at(j) > 0) {
                float deltaR1 = sqrt(pow(etagen->at(0)-rechiteta->at(j),2)+pow(phigen->at(0)-rechitphi->at(j),2));
                float deltaR2 = sqrt(pow(etagen->at(1)-rechiteta->at(j),2)+pow(phigen->at(1)-rechitphi->at(j),2));
                if(deltaR1 < 0.1 || deltaR2 < 0.1) rechitsum+=rechit->at(j);
            }
        }
        if (rechitsum > 1) h1->Fill(rechitsum);
    }
    h1->Draw();
    //c->SetLogy();
}
