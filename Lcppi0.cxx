//******* test of Reconstruction mctruth ****************
#include "Lcppi0Alg/Lcppi0.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IJobOptionsSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IMessageSvc.h"
#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/PropertyMgr.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EvtRecEvent/EvtRecVeeVertex.h"
#include "EvtRecEvent/EvtRecPi0.h"
#include "McTruth/McParticle.h"

#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/VertexFit.h"

#include "DstEvent/TofHitStatus.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"

#include "EvtRecEvent/EvtRecDTag.h"
#include "EvtRecEvent/EvtRecEtaToGG.h"
#include "MdcRecEvent/RecMdcKalTrack.h"

#include "TRandom.h"
#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"

#include "McDecayModeSvc/McDecayModeSvc.h"

#include "VertexFit/Helix.h"
#include "VertexFit/KinematicFit.h"
#include "ParticleID/ParticleID.h"
#include <vector>

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

typedef std::vector<int> Vint;
typedef std::vector<double> Vdou;
typedef std::vector<HepLorentzVector> Vp4;

using namespace Event;

const double xmass[4] = {0.000511, 0.105658, 0.1349766, 0.938272};//0=electron 1=muon 2=pi0 3=proton
const double mpi0 =  0.1349766;
const double mp = 0.9382723;


//*********************** Algorithm ***************************************
Lcppi0::Lcppi0(const std::string& name, ISvcLocator* pSvcLocator):Algorithm(name, pSvcLocator){
		m_ntot = 0;  m_debug = 0;
		declareProperty("Decay_ppi0", decay_ppi0   = 1);
		declareProperty("CheckTotal", m_checktotal = 1);
		declareProperty("HaveP",   m_HaveP   = 1);
		declareProperty("HavePi0",   m_HavePi0   = 1);
}


//********************** initialize ***************************************
StatusCode Lcppi0::initialize(){
		MsgStream log(msgSvc(), name());
		StatusCode scmgn = service ("MagneticFieldSvc",m_pIMF);
		if(scmgn!=StatusCode::SUCCESS) {
				log << MSG::ERROR << "Unable to open Magnetic field service"<<endreq;
		}
		log << MSG::INFO << "in initialize()" << endreq;
		StatusCode status;

cout << "before tree" << endl;

		if(decay_ppi0){
				NTuplePtr nt10(ntupleSvc(), "FILE1/ppi0");
				if(nt10) m_nt10 = nt10;

                else{
				  m_nt10 = ntupleSvc()->book("FILE1/ppi0", CLID_ColumnWiseTuple, "Lc+ Study");
				}
				if(m_nt10){

            status = m_nt10->addItem("n_Lc", n_Lc, 0, 1000);
            status = m_nt10->addItem("n_P", n_P, 0, 1000);
            status = m_nt10->addItem("n_Pi0", n_Pi0, 0, 1000);
            status = m_nt10->addItem("n_Gam", n_Gam, 0, 1000); 

            status = m_nt10->addIndexedItem("p4_Lc", n_Lc, 4, p4_Lc);
            status = m_nt10->addIndexedItem("p4_P", n_P, 4, p4_P);
            status = m_nt10->addIndexedItem("p4_Pi0", n_Pi0, 4, p4_Pi0);
            status = m_nt10->addIndexedItem("p4_Gam", n_Gam, 4, p4_Gam);

            status = m_nt10->addIndexedItem("Mbc", n_Lc, Mbc);
            status = m_nt10->addIndexedItem("deltaE", n_Lc, deltaE);
            status = m_nt10->addIndexedItem("chisquare", n_Pi0, chisquare);
            status = m_nt10->addIndexedItem("Mpi0", n_Pi0, Mpi0);            
            status = m_nt10->addIndexedItem("p_Vr", n_P, p_Vr);
            status = m_nt10->addIndexedItem("mdang", n_Gam, mdang);
            status = m_nt10->addIndexedItem("mpdang", n_Gam, mpdang);
        
        }
				else{
				  log << MSG::ERROR << "  Cannot book N-tuple10: " << long(m_nt10) << endreq;
                  return StatusCode::FAILURE;
				}
		}
		
		//Truth info write
/*      
		if(m_checktotal){
				NTuplePtr nt1(ntupleSvc(), "FILE1/tree_truth");
				if ( nt1 ) m_tuple1 = nt1;
				else {
						m_tuple1 = ntupleSvc()->book ("FILE1/tree_truth", CLID_ColumnWiseTuple, "tree_truth N-Tuple example");
             }
						if ( m_tuple1 )    {
								status = m_tuple1->addItem ("runNo",       m_runNo_);
								status = m_tuple1->addItem ("evtNo",       m_evtNo_);
								status = m_tuple1->addItem ("mode1",       m_mode1_);
								status = m_tuple1->addItem ("mode2",       m_mode2_);
								status = m_tuple1->addItem ("mode3",       m_mode3_);
								status = m_tuple1->addItem ("ndaughterAp",       m_ndaughterAp_, 0, 15);
								status = m_tuple1->addIndexedItem ("Ap_id", m_ndaughterAp_, m_Ap_id_);
								status = m_tuple1->addIndexedItem ("Ap_ptruth", m_ndaughterAp_,  4, m_Ap_ptruth_);
								status = m_tuple1->addItem ("ndaughterAm",       m_ndaughterAm_, 0, 15);
								status = m_tuple1->addIndexedItem ("Am_id", m_ndaughterAm_, m_Am_id_);
								status = m_tuple1->addIndexedItem ("Am_ptruth", m_ndaughterAm_,  4, m_Am_ptruth_);
						}
						else    { 
								log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) << endmsg;
								return StatusCode::FAILURE;
						}
				}*/
		log << MSG::INFO << "successfully return from initialize()" <<endmsg;
		return StatusCode::SUCCESS;
}
//End of initialize

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
StatusCode Lcppi0::execute(){
//    cout << "###############################" << endl;
		m_ntot = m_ntot + 1;
		MsgStream log(msgSvc(), name());
		log << MSG::INFO << "in execute()" << endreq;

		SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
		if (!eventHeader){
				log << MSG::FATAL << "Could not find Event Header" << endreq;
				return( StatusCode::FAILURE);
		}

		int runNo = eventHeader->runNumber();
		int eventNo = eventHeader->eventNumber();

		// ***** topology
		IMcDecayModeSvc* i_svc;
		StatusCode sc_DecayModeSvc = service("McDecayModeSvc", i_svc);
		if(sc_DecayModeSvc.isFailure()){
				log << MSG::FATAL << "Could not load McDecayModeSvc!" << endreq;
				return sc_DecayModeSvc;
		}
		m_svc = dynamic_cast<McDecayModeSvc*>(i_svc);

		//Truth info acquire
/*
		int runNo = eventHeader->runNumber();
		int eventNo = eventHeader->eventNumber();
		int mm_mode1=((eventHeader->flag1()/1000000))%1000;
		int mm_mode2=(eventHeader->flag1()/1000)%1000 ;
		int mm_mode3=eventHeader->flag1()%1000;


		// ***** topology
		IMcDecayModeSvc* i_svc;
		StatusCode sc_DecayModeSvc = service("McDecayModeSvc", i_svc);
		if(sc_DecayModeSvc.isFailure()){
				log << MSG::FATAL << "Could not load McDecayModeSvc!" << endreq;
				return sc_DecayModeSvc;
		}
		m_svc = dynamic_cast<McDecayModeSvc*>(i_svc);

		int M_pdgid[100];
		int M_motheridx[100];

		int M_pdgid_p[100];
		int M_motheridx_p[100];

		int M_pdgid_m[100];
		int M_motheridx_m[100];

		int numParticle = 0;
		int numParticle_p = 0;
		int numParticle_m = 0;

		int ndaughterAp=0; double	Ap_ptruth[15][4]; int Ap_id[15];
		for ( int aa = 0; aa < 15; aa++ ) 
				for ( int ll = 0; ll < 4; ll++ ) 
						Ap_ptruth[aa][ll]=0;
		for ( int aa = 0; aa < 15; aa++ ) 
				Ap_id[aa]=0;

		int	ndaughterAm=0; double	Am_ptruth[15][4]; int Am_id[15];
		for ( int aa = 0; aa < 15; aa++ ) 
				for ( int ll = 0; ll < 4; ll++ ) 
						Am_ptruth[aa][ll]=0;
		for ( int aa = 0; aa < 15; aa++ ) 
				Am_id[aa]=0;

     Vp4 gammar; gammar.clear();				
     Vp4 neutron; neutron.clear();
     Vp4 anti_neutron; anti_neutron.clear();
		if (eventHeader->runNumber()<0)
		{
				SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
				if (!mcParticleCol)
				{
						std::cout << "Could not retrieve McParticelCol" << std::endl;
						return StatusCode::FAILURE;
				}
				else {
						std::vector<int> pdgid;
						std::vector<int> motherindex;
						pdgid.clear();
						motherindex.clear();
//						bool Lmdc_P=false;       bool Lmdc_M=false;
						Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
            
            cout << "Truth info" << endl;
            cout << endl;
							for (iter_mc = mcParticleCol->begin(); iter_mc != mcParticleCol->end(); iter_mc++){
								if ((*iter_mc)->primaryParticle()) continue;
								if (!(*iter_mc)->decayFromGenerator()) continue;

								int pdg = (*iter_mc)->particleProperty();
								int motherpdg = ((*iter_mc)->mother()).particleProperty();
								int mmotherpdg = (((*iter_mc)->mother()).mother()).particleProperty();


           //Lambda_c+
        		if ((*iter_mc)->particleProperty()==4122){
				    int mode = m_svc->extract(*iter_mc, pdgid, motherindex);
    				numParticle_p = pdgid.size();
           // cout << "Decay list of  Lambda_c+" << endl;
  	       		for (int i=0; i!=  pdgid.size(); i++ ) {
					     	M_pdgid_p[i] = pdgid[i];
		  			 //  	if(!m_debug) cout<<"M_pdgid_p[i]="<<M_pdgid_p[i]<<endl;
						  	M_motheridx_p[i] = motherindex[i];
			  			 // if(!m_debug) cout<<"M_motheridx_p[i]="<<M_motheridx_p[i]<<endl;
                }  
					const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
					for(unsigned int ii = 0; ii < gc.size(); ii++) {
					if( gc[ii]->particleProperty() == -22) continue;
				    	Ap_id[ndaughterAp]=gc[ii]->particleProperty();
			  	for( int ll = 0; ll < 4; ll++ )
              Ap_ptruth[ndaughterAp][ll]=gc[ii]->initialFourMomentum()[ll];
						  ndaughterAp++;
						} 
          }

          //Pi0 from Lambda_c+
 					if(pdg==111&&motherpdg==4122){//pi0(pdg)->gam gam
					const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
							for(unsigned int ii = 0; ii < gc.size(); ii++) {
								if( gc[ii]->particleProperty() == -22) continue;
									Ap_id[ndaughterAp]=gc[ii]->particleProperty();
							for( int ll = 0; ll < 4; ll++ ) Ap_ptruth[ndaughterAp][ll]=gc[ii]->initialFourMomentum()[ll];
												ndaughterAp++;
										}// End of "gc.size() > 0" IF
								}
          //cout << endl; 
          //Lambda_c-
  					if ((*iter_mc)->particleProperty()==-4122){
	         // cout << "Decay list of  Lambda_c-" << endl;
  					int mode = m_svc->extract(*iter_mc, pdgid, motherindex);
						numParticle_m = pdgid.size();
								for (int i=0; i!=  pdgid.size(); i++ ) {
								M_pdgid_m[i] = pdgid[i];
							 if(!m_debug) cout<< M_pdgid_m[i]<<endl;
								M_motheridx_m[i] = motherindex[i];
								//if(!m_debug) cout<<"M_motheridx_m[i]="<<M_motheridx_m[i]<<endl;
										}
					const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
					for(unsigned int ii = 0; ii < gc.size(); ii++) {
					if( gc[ii]->particleProperty() == -22) continue;
				    	Am_id[ndaughterAm]=gc[ii]->particleProperty();
			  	for( int ll = 0; ll < 4; ll++ )
              Am_ptruth[ndaughterAm][ll]=gc[ii]->initialFourMomentum()[ll];
						  ndaughterAm++;
						} 
					}

            //Pi0 from Lambda_c-
         	if(pdg==111&&motherpdg==-4122){//pi0(pdg)->gam gam
					const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
							for(unsigned int ii = 0; ii < gc.size(); ii++) {
								if( gc[ii]->particleProperty() == -22) continue;
                if( gc[ii]->particleProperty() == 22)
                  { 
                    gammar.push_back(gc[ii]->initialFourMomentum());
                  }
									Am_id[ndaughterAm]=gc[ii]->particleProperty();
							for( int ll = 0; ll < 4; ll++ ) Am_ptruth[ndaughterAm][ll]=gc[ii]->initialFourMomentum()[ll];
												ndaughterAm++;
										}// End of "gc.size() > 0" IF
								}
        
  
         //Gammar from Lambda_c-
         if((pdg==22 && motherpdg != 111) || (pdg==22 && motherpdg == 111 && mmotherpdg != 4122))
         {
           gammar.push_back((*iter_mc)->initialFourMomentum());       
         }
         
         //neutron 
         if (pdg == 2112)
         {
           neutron.push_back((*iter_mc)->initialFourMomentum());
         }
         //anti-neutron 
         if (pdg == -2112)
         {
          anti_neutron.push_back((*iter_mc)->initialFourMomentum());
         }


         }
    }
    }
    
    cout << endl;

		if(m_checktotal){
				m_runNo_=runNo;
				m_evtNo_=eventNo;
				m_mode1_=mm_mode1;
				m_mode2_=mm_mode2;
				m_mode3_=mm_mode3;
				m_ndaughterAp_= ndaughterAp;
  
        double id[ndaughterAp];
        double p4[ndaughterAp][4];
				for ( int aa = 0; aa < ndaughterAp; aa++ )
          {        
             m_Ap_id_[aa]=Ap_id[aa];
             id[aa] = Ap_id[aa];
          }
				for ( int aa = 0; aa < ndaughterAp; aa++ ) 
          {
						for ( int ll = 0; ll < 4; ll++ ) 
							{
                	m_Ap_ptruth_[aa][ll]=Ap_ptruth[aa][ll];
                  p4[aa][ll] = Ap_ptruth[aa][ll];
              }
          }  
       
        cout << "p4 for all particles" << endl;
        cout << "4122" << ":\t("
             << p4[0][0] + p4[1][0] << ","
             << p4[0][1] + p4[1][1] << ","
             << p4[0][2] + p4[1][2] << ","
             << p4[0][3] + p4[1][3] << ")"
             << endl;
        for (int i = 0; i < ndaughterAp; i++)
          {
             cout << id[i] << ":\t("
                  << p4[i][0] << ","
                  << p4[i][1] << ","
                  << p4[i][2] << ","
                  << p4[i][3] << ")"
                  << endl;
          }
   
       //gammar from Lambda_c-
       cout << endl;
       cout << "p4 of Gammars from Lambda_c-" << endl;
       for (int i = 0; i < gammar.size(); i++)
          {
           cout << "22\t" << gammar[i] << endl;
          }
       cout << endl;
       //neutron
       cout << "p4 of neutron" << endl;
       for (int i = 0; i < neutron.size(); i++)
          {
           cout << "2112\t" << neutron[i] << endl;
          }
       cout << endl;
      //anti-neutron
       cout << "p4 of anti-neutron" << endl;
       for (int i = 0; i < anti_neutron.size(); i++)
          {
           cout << "-2112\t" << anti_neutron[i] << endl;
          }
       cout << endl;

        	HepLorentzVector p4_p;
        	HepLorentzVector  p4_pi0;
        	HepLorentzVector  p4_Lc;

        for ( int jj = 0; jj < 4; jj++)
           {
             p4_p[jj] = Ap_ptruth[0][jj];
             p4_pi0[jj] = Ap_ptruth[1][jj];
             p4_Lc[jj] = Ap_ptruth[0][jj] + Ap_ptruth[1][jj];
           }

        cout << endl;
        cout << "Lambda_c+_p4_truth"<< p4_Lc << endl;
        cout << "Proton p4_truth: " << p4_p << endl;
        cout << "  Pi0_p4_truth:  " << p4_pi0 << endl;        

				m_ndaughterAm_= ndaughterAm;
				for ( int aa = 0; aa < ndaughterAm; aa++ ) m_Am_id_[aa]=Am_id[aa];
				for ( int aa = 0; aa < ndaughterAm; aa++ ) 
						for ( int ll = 0; ll < 4; ll++ ) 
								m_Am_ptruth_[aa][ll]=Am_ptruth[aa][ll];

				m_tuple1->write();
		}

      cout << "End of Truth info" << endl;
      cout << endl;*/


		// >>>>>>>>>>>>>>>>>>>>  Selection and Reconstruction begin from here <<<<<<<<<<<<<<<<<<<<<<    //before this the tracks have been found and rebuilded
		                                                                                                //related properties have been stored.
		SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);   //a class like Vector, store if the correspoing track is charged
		SmartDataPtr<EvtRecTrackCol> evtRecTrackCol( eventSvc(), EventModel::EvtRec::EvtRecTrackCol);  //a class like Vector, store the details of every track

		Hep3Vector xorigin(0,0,0);
		IVertexDbSvc* vtxsvc;
		Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);   //the coordinates of vertex of event?
		if(vtxsvc->isVertexValid()){
				double* dbv = vtxsvc->PrimaryVertex(); //coordinates of primary vertex
			    double*  vv = vtxsvc->SigmaPrimaryVertex();  
				xorigin.setX(dbv[0]);
				xorigin.setY(dbv[1]);
				xorigin.setZ(dbv[2]);
		}


	Vp4 p4pp;  p4pp.clear();
    Vdou pVr; pVr.clear();
		for(int i = 0; i < evtRecEvent->totalCharged(); i++){
				EvtRecTrackIterator itTrk = evtRecTrackCol->begin() + i;
				RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack(); //store information of particle like charge
				RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//relatted to pid and calculation of properties of particles

				//judgement of the charged-Track
				if(!isGoodTrack(*itTrk)) continue;

				//Particle IDentify for each charged-Track
				preparePID(*itTrk);
				if(mdcTrk->charge()==1){
                  double Vr = 10.;
				  if(m_HaveP&&isproton(*itTrk, Vr)){
			    	mdcKalTrk->setPidType(RecMdcKalTrack::proton);
					p4pp.push_back(mdcKalTrk->p4(xmass[3]));  //with the information of mass of proton and from track, calculate the p4 of proton.i
                    pVr.push_back(Vr);
				  }
				}
		}
     int nP = p4pp.size();
     if(nP < 1) return StatusCode::SUCCESS;

		// *** select good shower *** //  
	Vint iGam; iGam.clear();
	Vp4 p4Gam; p4Gam.clear();
    Vdou Dang; Dang.clear();
    Vdou mpDang; mpDang.clear();
		for(int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++){   //tracks are divided into two parts, one is charged, the other is not
				EvtRecTrackIterator itTrk = evtRecTrackCol->begin() + i;                     //when select protons, we seek it in the charged part, now we seek for gammer
        double dang = 200.;
        double mpdang = 200.;
				if(!isGoodShower(*itTrk, dang, mpdang)) continue;                          //which is not charged, start from totalCharged, end at totalTracks.
		     iGam.push_back(i);
				RecEmcShower* emcTrk = (*itTrk)->emcShower();                                //evtRecTrackCol is some sort of array recording all tracks. from i=0 to
				p4Gam.push_back(getP4(emcTrk,xorigin)); 
        Dang.push_back(dang);
        mpDang.push_back(mpdang);
        //why we need xorigin to calculate p4 of gammer                              //i = totalCharged() it stores charged tracks, from i = totalCharged() to 
		}   //it is the coordinate of primary vertex, if conjecture if right should we use the coordinates of vertex where pi0->2gam
                                                                                     //i = totalTracks() it stores uncharged tracks which are recored by EMC.
		int nGam = p4Gam.size();
		if(nGam<1) return StatusCode::SUCCESS;

		Vdou massPi0; massPi0.clear();
		Vdou chisPi0; chisPi0.clear();
		Vp4 p4Pi0; p4Pi0.clear();
		Vp4 p4Pi01c; p4Pi01c.clear();
		int nPi0;
    
		if(m_HavePi0){
				// Loop each gamma pair, check if it is a pi0
				for(int i = 0; i < nGam-1; i++) 
				{
						for(int j = i+1; j < nGam; j++) 
						{
								EvtRecTrackIterator itTrki = evtRecTrackCol->begin() + iGam[i];
								RecEmcShower* shr1 = (*itTrki)->emcShower();

								EvtRecTrackIterator itTrkj = evtRecTrackCol->begin() + iGam[j];
								RecEmcShower* shr2 = (*itTrkj)->emcShower();

								HepLorentzVector p4_pi0(0,0,0,0),p4_pi0_1c(0,0,0,0);
								double pi0_mass;
								double pi0_chis;
								if(isGoodpi0(shr1,shr2,pi0_mass,p4_pi0,pi0_chis,p4_pi0_1c)){
										massPi0.push_back(pi0_mass);
										chisPi0.push_back(pi0_chis);
										p4Pi0.push_back(p4_pi0);
										p4Pi01c.push_back(p4_pi0_1c);
              					}

						}
				}
				nPi0 = p4Pi01c.size();
				if(nPi0 < 1) return StatusCode::SUCCESS;
		}


		//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
		double ECMS = 2.3;
		if(runNo <= 0) ECMS = 2.3;
		if(runNo >=43716 && runNo<= 47185) ECMS = 2.3;


Vp4 p4Lc; p4Lc.clear();
Vdou delta_E; delta_E.clear();
Vdou mbc; mbc.clear();
		for(int k = 0; k < nP; k++){
			for(int j = 0; j < nPi0; j++){      
                    HepLorentzVector  m_p4Lc  =  p4pp[k] + p4Pi01c[j];
                    m_p4Lc.boost(-0.011,0,0);
					double mbc2 = ECMS*ECMS -m_p4Lc.v().mag2();
					double m_Lc = mbc2 > 0 ? sqrt(mbc2) : -10;
                    p4Lc.push_back(m_p4Lc);
                    mbc.push_back(m_Lc);
	               	delta_E.push_back(m_p4Lc.t() - ECMS);
				}
              
		}
      int nLc = p4Lc.size();
      if(nLc < 1) return StatusCode::SUCCESS;
				

//store variables
//numbers
       n_Lc = nLc;
       n_P = nP;
       n_Pi0 = nPi0;
       n_Gam = nGam;
//p4
       for(int i = 0; i < 4; i++)
         {  
           for(int j = 0; j < nLc; j++)
            {
               p4_Lc[j][i] = p4Lc[j][i];    
            }

           for(int j = 0; j < nP; j++)
            {
               p4_P[j][i] = p4pp[j][i];    
            }

           for(int j = 0; j < nPi0; j++)
            {
               p4_Pi0[j][i] = p4Pi01c[j][i];    
            }

           for(int j = 0; j < nGam; j++)
            {
               p4_Gam[j][i] = p4Gam[j][i];    
            }
         }
//other
       for(int i = 0; i < nLc; i++)
        {
           Mbc[i] = mbc[i];
           deltaE[i] = delta_E[i];
        }
       
       for(int i = 0; i < nPi0; i++)
        {
           chisquare[i] = chisPi0[i];
           Mpi0[i] = massPi0[i]; 
        }
         
       for(int i = 0; i < nP; i++)
        {
          p_Vr[i] = pVr[i];
        }

       for(int i = 0; i < nGam; i++)
        {
          mdang[i] = Dang[i];
          mpdang[i] = mpDang[i];
        } 
//write
        m_nt10->write();

		return StatusCode::SUCCESS;
}


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

StatusCode Lcppi0::finalize() {
		MsgStream log(msgSvc(), name());
		log << MSG::INFO << "in finalize()" << endreq;
		return StatusCode::SUCCESS;
}


// * * * * * * * * * * * * * * * * * * * Sub Program * * * * * * * * * * * * * * * * * * * //
HepLorentzVector Lcppi0::getP4(RecEmcShower* gTrk, Hep3Vector origin){
		Hep3Vector Gm_Vec(gTrk->x(), gTrk->y(), gTrk->z());
		Hep3Vector Gm_Mom = Gm_Vec - origin;
		Gm_Mom.setMag(gTrk->energy());
		HepLorentzVector pGm(Gm_Mom, gTrk->energy());
		return pGm;
}

bool Lcppi0::isGoodTrack(EvtRecTrack* trk){
		double m_vz0cut = 10.0;
		double m_vr0cut = 1.0;
		double m_CosThetaCut = 0.93;

		if(!trk->isMdcKalTrackValid()){
				return false;
		}
		Hep3Vector xorigin(0,0,0);
		IVertexDbSvc*  vtxsvc;
		Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
		if(vtxsvc->isVertexValid()){
				double* dbv = vtxsvc->PrimaryVertex();
				double* vv = vtxsvc->SigmaPrimaryVertex();
				xorigin.setX(dbv[0]);
				xorigin.setY(dbv[1]);
				xorigin.setZ(dbv[2]);
		}
		RecMdcTrack *mdcTrk = trk->mdcTrack();
		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.); // the initial point for MDC recosntruction
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
		VFHelix helixip(point0,a,Ea);
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
		double Rvz0=vecipa[3];         //the nearest distance to IP in z direction
		double costheta=cos(mdcTrk->theta());
		if(fabs(Rvz0) < m_vz0cut && fabs(Rvxy0)< m_vr0cut && fabs(costheta)<m_CosThetaCut) return true;
		return false;
}

void Lcppi0::preparePID(EvtRecTrack* track){
		//PID for proton
		ParticleID *pid = ParticleID::instance();
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(8);
		pid->setRecTrack(track);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTof()); // use PID sub-system
		//    pid->identify(pid->onlyPion() | pid->onlyKaon());
		pid->identify(pid->onlyPion()| pid->onlyKaon() | pid->onlyProton());//change
		pid->calculate();
		m_prob[0]=pid->probElectron();
		m_prob[1]=pid->probMuon();
		m_prob[2]=pid->probPion();
		m_prob[3]=pid->probKaon();
		m_prob[4]=pid->probProton();
		if(!(pid->IsPidInfoValid())){
				for(int i=0; i<5; i++) m_prob[i]=-99;
		}
}


bool Lcppi0::isproton(EvtRecTrack* trk, double &Vr){
		if (m_prob[4]>=0.00 && m_prob[4]>=m_prob[3] && m_prob[4]>=m_prob[2])
     { 
   		Hep3Vector xorigin(0,0,0);
	  	IVertexDbSvc*  vtxsvc;
		  Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
  		if(vtxsvc->isVertexValid()){
	  			double* dbv = vtxsvc->PrimaryVertex();
		  		double* vv = vtxsvc->SigmaPrimaryVertex();
			  	xorigin.setX(dbv[0]);
				  xorigin.setY(dbv[1]);
  				xorigin.setZ(dbv[2]);
	  	}
		RecMdcTrack *mdcTrk = trk->mdcTrack();
  		HepVector a = mdcTrk->helix();
	  	HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.); // the initial point for MDC recosntruction
  		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
	  	VFHelix helixip(point0,a,Ea);
		helixip.pivot(IP);
  		HepVector vecipa = helixip.a();
    	Vr = fabs(vecipa[0]);  //the nearest distance to IP in xy plane 
        return true;
     }
		return false;
}

bool Lcppi0::isGoodShower(EvtRecTrack* trk, double &dang, double &mpdang){
		double m_maxCosThetaBarrel = 0.80;
		double m_minCosThetaEndcap = 0.86;
		double m_maxCosThetaEndcap = 0.92;
		double m_minBarrelEnergy = 0.025;
		double m_minEndcapEnergy = 0.050;
		double m_gammthCut = 8.0; //maximum angle between showers and charged track
        double m_pCut = 20.;

		if(!trk->isEmcShowerValid()){
				return false;
		}
		Hep3Vector xorigin(0,0,0);
		IVertexDbSvc*  vtxsvc;
		Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
		if(vtxsvc->isVertexValid()){
				double* dbv = vtxsvc->PrimaryVertex();
				double*  vv = vtxsvc->SigmaPrimaryVertex();
				xorigin.setX(dbv[0]);
				xorigin.setY(dbv[1]);
				xorigin.setZ(dbv[2]);
		}
		SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
		SmartDataPtr<EvtRecTrackCol> evtRecTrackCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
		RecEmcShower *emcTrk = trk->emcShower();
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());   //position of shower
		HepLorentzVector shP4 = getP4(emcTrk,xorigin);
		double cosThetaSh = shP4.vect().cosTheta();
		double eraw = emcTrk->energy();  
		double getTime = emcTrk->time(); //in unit of 50ns
		if(getTime>14||getTime<0) return false; //cut time
		if(!((fabs(cosThetaSh)<m_maxCosThetaBarrel && eraw>m_minBarrelEnergy)||((fabs(cosThetaSh)>m_minCosThetaEndcap) && (fabs(cosThetaSh)<m_maxCosThetaEndcap) && (eraw>m_minEndcapEnergy))))  return false;
	    dang = 200.;
        mpdang = 200.;
  	for(int j = 0; j < evtRecEvent->totalCharged(); j++){
				EvtRecTrackIterator jtTrk = evtRecTrackCol->begin() + j; 
        RecMdcTrack *mdcTrk = (*jtTrk)->mdcTrack();  //track in mdc;
				if(!(*jtTrk)->isExtTrackValid()) continue;   
				RecExtTrack *extTrk = (*jtTrk)->extTrack();
				if(extTrk->emcVolumeNumber() == -1) continue;
				Hep3Vector extpos = extTrk->emcPosition();   //position of 
				double angd = extpos.angle(emcpos);
        preparePID(*jtTrk);
        double Vr = 10.;
        if(mdcTrk->charge() == -1 && isproton(*jtTrk, Vr))
        mpdang = (mpdang > angd) ? angd : mpdang;
				if(angd < dang) dang = angd;
		}
		if(dang>=200) return false;               //if there is no ExtTrack or no emc, still return false?
		dang = dang * 180 / (CLHEP::pi);
        mpdang = mpdang * 180 / (CLHEP::pi);
		return true;
}

bool Lcppi0::isGoodpi0(RecEmcShower *shr1,RecEmcShower *shr2,double& pi0_mass,HepLorentzVector& p4_pi0,double& pi0_chis,HepLorentzVector& p4_pi0_1c){

		Hep3Vector xorigin(0,0,0);
		IVertexDbSvc*  vtxsvc;
		Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
		if(vtxsvc->isVertexValid()){
				double* dbv = vtxsvc->PrimaryVertex();
				double*  vv = vtxsvc->SigmaPrimaryVertex();
				xorigin.setX(dbv[0]);
				xorigin.setY(dbv[1]);
				xorigin.setZ(dbv[2]);
		}

		pi0_mass=-100;
		pi0_chis=-100;

		HepLorentzVector g1P4 = getP4(shr1,xorigin);
		HepLorentzVector g2P4 = getP4(shr2,xorigin);
		p4_pi0 = g1P4 + g2P4;
		pi0_mass = p4_pi0.m();

		double xmpi0=0.134976;

		KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
		kmfit->init();
		kmfit->setIterNumber(5);
		kmfit->AddTrack(0, 0.0, shr1);
		kmfit->AddTrack(1, 0.0, shr2);
		kmfit->AddResonance(0, xmpi0, 0, 1);

	//	if(!olmdq) return false;
    if(kmfit->Fit(0))
	  {
    	kmfit->BuildVirtualParticle(0);
		  pi0_chis = kmfit->chisq(0);
      }
    else
    pi0_chis = 100.; //if fit falls, let chis = 100;

	p4_pi0_1c=kmfit->pfit(0)+kmfit->pfit(1);
	return true;
}


//******************************************************************
