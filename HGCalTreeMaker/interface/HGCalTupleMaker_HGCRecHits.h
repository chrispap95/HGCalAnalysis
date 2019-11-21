#ifndef HGCalTupleMaker_HGCRecHits_h
#define HGCalTupleMaker_HGCRecHits_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/transform.h"

// HGCAL objects
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

// Geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"

// Used:
// http://cmslxr.fnal.gov/source/AnalysisAlgos/SiStripClusterInfoProducer/plugins/SiStripProcessedRawDigiProducer.cc
// as an example

//-----------------------------------------------------
// Class definition
//-----------------------------------------------------
class HGCalTupleMaker_HGCRecHits : public edm::EDProducer {
protected:
    std::vector<edm::InputTag> m_HGCRecHitsTags;
    std::vector<edm::EDGetTokenT<HGCRecHitCollection>> m_HGCRecHitsTokens;
    std::vector<std::string> m_geometrySource;

    const std::string     m_prefix;
    const std::string     m_suffix;

    bool debug=false;

    //HGC Geometry
    std::vector<const HGCalDDDConstants*> hgcCons_;
    std::vector<const HGCalGeometry*>     hgcGeometry_;
    const HcalDDDSimConstants*            hcCons_;
    const HcalDDDRecConstants*            hcConr_;
    const CaloSubdetectorGeometry*        hcGeometry_;

    edm::ESHandle<CaloGeometry> geometry;
    std::map<uint32_t, HepGeom::Transform3D> transMap_;

    void produce( edm::Event & iEvent, const edm::EventSetup & iSetup) {
        //-----------------------------------------------------
        // Prepare to put things into event
        //-----------------------------------------------------
        loadAlgo();

        //-----------------------------------------------------
        // edm::Handles
        //-----------------------------------------------------
        edm::Handle<HGCRecHitCollection> HGCRecHits;

        //-----------------------------------------------------
        // Loop over HGCRecHitCollections
        //-----------------------------------------------------
        for( typename std::vector<edm::EDGetTokenT<HGCRecHitCollection> >::const_iterator
          token = m_HGCRecHitsTokens.begin(); token != m_HGCRecHitsTokens.end(); ++token ) {
            unsigned index(token - m_HGCRecHitsTokens.begin());
            HGCRecHits.clear();

            //---------------------------------------------------
            // Geometry & looping over rechits
            //---------------------------------------------------
            std::string nameDetector_ = m_geometrySource[index];
            if (debug) std::cout << nameDetector_ << std::endl;

            //-----------------------------------------------------
            // getByToken
            //-----------------------------------------------------
            iEvent.getByToken(*token, HGCRecHits);
            if( HGCRecHits.isValid() && !HGCRecHits->empty() ) {
                edm::LogInfo("Input found") << m_HGCRecHitsTags.at(index);
            } else {
                edm::LogInfo("Input not found") << m_HGCRecHitsTags.at(index);
                return;
            }

//******************************************************************************
            //* HCal's geometry
            //* int cell, type, sector, subsector, layer, zside;
            //* int type, zside;
            //* int cellU=-100, cellV=-100, waferU=-100, waferV=-100;
            //* int ieta=-100, iphi=-100, ietaAbs=-100;
            //* int subdet(0);
            //* HepGeom::Point3D<float> gcoord;

            if (nameDetector_ == "HCal") {
                edm::ESHandle<CaloGeometry> geom;
                iSetup.get<CaloGeometryRecord>().get(geom);
                if (!geom.isValid()) {
                    edm::LogWarning("HGCalTupleMaker_HGCRecHits")
                    << "Cannot get valid HGCalGeometry Object for " << nameDetector_;
                    return;
                }
                const CaloGeometry* geom0 = geom.product();

                for (const auto & it : *(HGCRecHits.product())) {
                    DetId detId = it.id();
                    //* subdet = detId.subdetId();
                    int ilayer = HcalDetId(detId).depth();
                    //* auto cellGeometry = geometry->getSubdetectorGeometry(detId)->getGeometry(detId);
                    //* gcoord = HepGeom::Point3D<float>(cellGeometry->getPosition().x(),
                    //* cellGeometry->getPosition().y(),
                    //* cellGeometry->getPosition().z());
                    run(detId, ilayer, index, geom0, &it);
                }
            } else {
                // HGCal's geometry
                edm::ESHandle<HGCalGeometry> geom;
                iSetup.get<IdealGeometryRecord>().get(nameDetector_, geom);
                if (!geom.isValid()) {
                    edm::LogWarning("HGCalTupleMaker_HGCRecHits")
                    << "Cannot get valid HGCalGeometry Object for " << nameDetector_;
                    return;
                }
                const HGCalGeometry* geom0 = geom.product();
                HGCalGeometryMode::GeometryMode mode = geom0->topology().geomMode();
                int geomType = (((mode == HGCalGeometryMode::Hexagon8) ||
                (mode == HGCalGeometryMode::Hexagon8Full)) ? 1 :
                ((mode == HGCalGeometryMode::Trapezoid) ? 2 : 0));

                for (const auto & it : *(HGCRecHits.product())) {
                    DetId detId = it.id();
                    //KH int ilayer   = HGCalDetId(detId).layer();
                    int ilayer  = ((geomType == 0) ? HGCalDetId(detId).layer() :
                        ((geomType == 1) ? HGCSiliconDetId(detId).layer() :
                        HGCScintillatorDetId(detId).layer()));
                    //* auto cellGeometry = geometry->getSubdetectorGeometry(detId)->getGeometry(detId);
                    //*     gcoord = HepGeom::Point3D<float>(cellGeometry->getPosition().x(),
                    //*     cellGeometry->getPosition().y(),
                    //*     cellGeometry->getPosition().z());

                    //* HGCSiliconDetId detIdSi = HGCSiliconDetId(detId);
                    //* subdet           = ForwardEmpty;
                    //* cellU            = detIdSi.cellU();
                    //* cellV            = detIdSi.cellV();
                    //* waferU           = detIdSi.waferU();
                    //* waferV           = detIdSi.waferV();
                    //* ieta             = detIdSi.ietaAbs();
                    //* iphi             = detIdSi.iphi();
                    //* type             = detIdSi.type();
                    //* layer            = detIdSi.layer();
                    //* zside            = detIdSi.zside();
                    //* v_ieta   -> push_back ( ieta   );
                    //* v_iphi   -> push_back ( iphi   );
                    //* v_cellu  -> push_back ( cellU  );
                    //* v_cellv  -> push_back ( cellV  );
                    //* v_waferu -> push_back ( waferU );
                    //* v_waferv -> push_back ( waferV );

                    run(detId, ilayer, index, geom0, &it);
                }
            } // Hcal or HGcal for geometry
        } // Looping over different rechit collections

        //-----------------------------------------------------
        // Put things into the event
        //-----------------------------------------------------
        dumpAlgo(iEvent);
    }

public:
    HGCalTupleMaker_HGCRecHits(const edm::ParameterSet& iConfig) :
      m_HGCRecHitsTags  (iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("source")),
      m_HGCRecHitsTokens(edm::vector_transform(m_HGCRecHitsTags, [this](edm::InputTag const & tag){
        return consumes<HGCRecHitCollection>(tag);})),
      m_geometrySource  (iConfig.getUntrackedParameter<std::vector<std::string> >("geometrySource")),
      m_prefix          (iConfig.getUntrackedParameter<std::string>  ("Prefix")),
      m_suffix          (iConfig.getUntrackedParameter<std::string>  ("Suffix")) {
        produces<std::vector<float> > ( m_prefix + "Energy" + m_suffix );
        //* produces<std::vector<float> > ( m_prefix + "Time"   + m_suffix );
        //* produces<std::vector<int  > > ( m_prefix + "Subdet" + m_suffix );
        produces<std::vector<int  > > ( m_prefix + "Layer"  + m_suffix );
        produces<std::vector<int  > > ( m_prefix + "Index"  + m_suffix );
        produces<std::vector<float> > ( m_prefix + "Eta"    + m_suffix );
        produces<std::vector<float> > ( m_prefix + "Phi"    + m_suffix );
        //* produces<std::vector<int  > > ( m_prefix + "IEta"   + m_suffix );
        //* produces<std::vector<int  > > ( m_prefix + "IPhi"   + m_suffix );
        //* produces<std::vector<int  > > ( m_prefix + "CellU"  + m_suffix );
        //* produces<std::vector<int  > > ( m_prefix + "CellV"  + m_suffix );
        //* produces<std::vector<int  > > ( m_prefix + "WaferU" + m_suffix );
        //* produces<std::vector<int  > > ( m_prefix + "WaferV" + m_suffix );
        produces<std::vector<float> > ( m_prefix + "Posx"   + m_suffix );
        produces<std::vector<float> > ( m_prefix + "Posy"   + m_suffix );
        produces<std::vector<float> > ( m_prefix + "Posz"   + m_suffix );
    }

    std::unique_ptr<std::vector<float> > eta;
    std::unique_ptr<std::vector<float> > phi;
    std::unique_ptr<std::vector<int  > > layer;
    std::unique_ptr<std::vector<float> > energy;
    std::unique_ptr<std::vector<float> > posx;
    std::unique_ptr<std::vector<float> > posy;
    std::unique_ptr<std::vector<float> > posz;
    std::unique_ptr<std::vector<int  > > v_index;
    //* std::unique_ptr<std::vector<int  > > v_subdet;

    //* std::unique_ptr<std::vector<int  > > v_cell; //
    //* std::unique_ptr<std::vector<int  > > v_sector; // wafer
    //* std::unique_ptr<std::vector<int  > > v_subsector; // type
    //* std::unique_ptr<std::vector<int  > > v_layer;
    //* std::unique_ptr<std::vector<int  > > v_zside;

    //* std::unique_ptr<std::vector<float> > v_posx;
    //* std::unique_ptr<std::vector<float> > v_posy;
    //* std::unique_ptr<std::vector<float> > v_posz;
    //* std::unique_ptr<std::vector<float> > v_eta;
    //* std::unique_ptr<std::vector<float> > v_phi;

    //* std::unique_ptr<std::vector<int  > > v_ieta;
    //* std::unique_ptr<std::vector<int  > > v_iphi;
    std::unique_ptr<std::vector<int  > > v_cellu;
    std::unique_ptr<std::vector<int  > > v_cellv;
    std::unique_ptr<std::vector<int  > > v_waferu;
    std::unique_ptr<std::vector<int  > > v_waferv;

protected:
    void loadAlgo(){
        eta      = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
        phi      = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
        layer    = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
        energy   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
        posx     = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
        posy     = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
        posz     = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
        v_index  = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
        //* v_ieta   = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
        //* v_iphi   = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
        v_cellu  = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
        v_cellv  = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
        v_waferu = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
        v_waferv = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());*/
    }

    void dumpAlgo( edm::Event & iEvent ){
        iEvent.put( move(eta      ), m_prefix + "Eta"    + m_suffix );
        iEvent.put( move(phi      ), m_prefix + "Phi"    + m_suffix );
        iEvent.put( move(layer    ), m_prefix + "Layer"  + m_suffix );
        iEvent.put( move(energy   ), m_prefix + "Energy" + m_suffix );
        iEvent.put( move(posx     ), m_prefix + "Posx"   + m_suffix );
        iEvent.put( move(posy     ), m_prefix + "Posy"   + m_suffix );
        iEvent.put( move(posz     ), m_prefix + "Posz"   + m_suffix );
        iEvent.put( move(v_index  ), m_prefix + "Index"  + m_suffix );
        //* iEvent.put( move(v_ieta   ), m_prefix + "IEta"   + m_suffix );
        //* iEvent.put( move(v_iphi   ), m_prefix + "IPhi"   + m_suffix );
        iEvent.put( move(v_cellu  ), m_prefix + "CellU"  + m_suffix );
        iEvent.put( move(v_cellv  ), m_prefix + "CellV"  + m_suffix );
        iEvent.put( move(v_waferu ), m_prefix + "WaferU" + m_suffix );
        iEvent.put( move(v_waferv ), m_prefix + "WaferV" + m_suffix );
    }

    template<class T1, class T2>
    void run(DetId & detId, int ilayer, int index, const T1* geom, T2 it) {
        GlobalPoint global = geom->getPosition(detId);

        if (debug) std::cout << ilayer << " " << global.z() << std::endl;

        eta      -> push_back ( global.eta()     );
        phi      -> push_back ( global.phi()     );
        layer    -> push_back ( ilayer           );
        energy   -> push_back ( it->energy()     );
        posx     -> push_back ( global.x()       );
        posy     -> push_back ( global.y()       );
        posz     -> push_back ( global.z()       );
        v_index  -> push_back ( index            );
        //* v_ieta   -> push_back ( global.ietaAbs() );
        //* v_iphi   -> push_back ( global.iphi()    );
        //* v_cellu  -> push_back ( global.cellU()   );
        //* v_cellv  -> push_back ( global.cellV()   );
        //* v_waferu -> push_back ( global.waferU()  );
        //* v_waferv -> push_back ( global.waferV()  );*/
    }
};

#endif
