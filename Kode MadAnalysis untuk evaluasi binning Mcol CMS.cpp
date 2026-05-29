#include "SampleAnalyzer/User/Analyzer/mcoltaumuCMS.h"
using namespace MA5;
using namespace std;
#include "TLorentzVector.h"

double binmue0j[22] = {40,80,100,125,150,175,200,235,270,315,350,400,450,500,560,630,700,780,860,950,1040,1200};
int nmue0j[22] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int nevent = 0;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool mcoltaumuCMS::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  cout << "BEGIN Initialization" << endl;
  // initialize variables, histos
  cout << "END   Initialization" << endl;
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void mcoltaumuCMS::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;
  // saving histos
  cout << "mue0jdist = {"; 
  for(int i=0; i < 21; i++)
  {
	cout << "{" << binmue0j[i] << ", " << nmue0j[i] << "}, ";
  }
  cout << "{" << binmue0j[21] << ", " << nmue0j[21] << "}}" << endl;
  
  ofstream myfile;
  myfile.open ("output.txt");
  for(int i=0; i < 22; i++)
  {
	myfile << nmue0j[i] << endl;
  }
  myfile.close();
  cout << "END   Finalization" << endl;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool mcoltaumuCMS::Execute(SampleFormat& sample, const EventFormat& event)
{
	// Event weight
  double myWeight=1.;
  if (!Configuration().IsNoEventWeight() && event.mc()!=0) myWeight=event.mc()->weight();

  Manager()->InitializeForNewEvent(myWeight);

	vector<TLorentzVector> electrons;
	vector<TLorentzVector> muons;
	vector<int> eleccharge;
	vector<int> mucharge;
	int ne = 0;
	int nmu = 0;
	int ntau = 0;
	int ne5 = 0;
	int nmu10 = 0;
	int ntau20 = 0;
	int njet = 0;
	
  if (event.rec()!=0)
  {
	  //Electron
	  for (MAuint32 i=0;i<event.rec()->electrons().size();i++)
	{
      const RecLeptonFormat& elec = event.rec()->electrons()[i];
      if(elec.pt() > 5 && abs(elec.eta()) < 2.5)
	  {
		  ne5++;
	  }
	  if(elec.pt() > 10 && abs(elec.eta()) < 2.4)
	  {
		  ne++;
		  TLorentzVector temp;
		  temp.SetPxPyPzE(elec.px(),elec.py(),elec.pz(),elec.e());
		  electrons.push_back(temp);
		  eleccharge.push_back(elec.charge());
	  }  
    }
	
	// Looking through the reconstructed muon collection
    for (MAuint32 i=0;i<event.rec()->muons().size();i++)
    {
      const RecLeptonFormat& mu = event.rec()->muons()[i];
      if(mu.pt() > 10 && abs(mu.eta()) < 2.5)
	  {
		  nmu10++;
	  }
	  if(mu.pt() > 60 && abs(mu.eta()) < 2.4)
	  {
		  nmu++;
		  TLorentzVector temp;
		  temp.SetPxPyPzE(mu.px(),mu.py(),mu.pz(),mu.e());
		  muons.push_back(temp);
		  mucharge.push_back(mu.charge());
	  } 
    }
	
	// Looking through the reconstructed hadronic tau collection
    for (MAuint32 i=0;i<event.rec()->taus().size();i++)
    {
      const RecTauFormat& tau = event.rec()->taus()[i];
      if(tau.pt() > 20 && abs(tau.eta()) < 2.5)
	  {
		  ntau20++;
	  }
    }
	
	// Looking through the reconstructed jet collection
    for (MAuint32 i=0;i<event.rec()->jets().size();i++)
    {
      const RecJetFormat& jet = event.rec()->jets()[i];
      if(jet.pt() > 30 && abs(jet.eta()) < 4.7)
	  {
		  njet++;
	  }
    }
	
	// Transverse missing energy (MET)
    TLorentzVector met;
	met.SetXYZM	(event.rec()->MET().px(),event.rec()->MET().py(),0,0);
	
	if(ne == 1 && nmu == 1 && ne5 < 2 && nmu10 < 2 && ntau20 < 1 && njet == 1)
	{
		if(electrons[0].DeltaR(muons[0]) > 0.3 && abs(electrons[0].DeltaPhi(met)) < 0.7 && abs(electrons[0].DeltaPhi(muons[0])) > 2.2)
		{
			TLorentzVector tavis = electrons[0];
			double mvis = (tavis+muons[0]).M();
			double ptnuest = (tavis.Px()*met.Px()+tavis.Py()*met.Py())/tavis.Pt();
			double xtvis = tavis.Pt()/(tavis.Pt()+ptnuest);
			double mcol = mvis/sqrt(xtvis);
	
			if(mcol > binmue0j[0] && mcol < binmue0j[21])
			{
				for(int i; i < 21; i++)
				{
					if(mcol > binmue0j[i] && mcol < binmue0j[i+1])
					{
						nmue0j[i]++;
					}
				}
			}
			else if(mcol > binmue0j[21])
			{
				nmue0j[21]++;
			}
		}
	}
	
	
	
	
  }
  
  
  // ***************************************************************************
  // Example of analysis with reconstructed objects
  // Concerned samples : 
  //   - LHCO samples
  //   - LHE/STDHEP/HEPMC samples after applying jet-clustering algorithm
  // ***************************************************************************
  /*
  // Event weight
  double myWeight=1.;
  if (!Configuration().IsNoEventWeight() && event.mc()!=0) myWeight=event.mc()->weight();

  Manager()->InitializeForNewEvent(myWeight);

  if (event.rec()!=0)
  {
    cout << "---------------NEW EVENT-------------------" << endl;

    // Looking through the reconstructed electron collection
    for (MAuint32 i=0;i<event.rec()->electrons().size();i++)
    {
      const RecLeptonFormat& elec = event.rec()->electrons()[i];
      cout << "----------------------------------" << endl;
      cout << "Electron" << endl;
      cout << "----------------------------------" << endl;
      cout << "index=" << i+1 
                << " charge=" << elec.charge() << endl;
      cout << "px=" << elec.px()
                << " py=" << elec.py()
                << " pz=" << elec.pz()
                << " e="  << elec.e()
                << " m="  << elec.m() << endl;
      cout << "pt=" << elec.pt() 
                << " eta=" << elec.eta() 
                << " phi=" << elec.phi() << endl;
      cout << "pointer address to the matching MC particle: " 
                << elec.mc() << endl;
      cout << endl;
    }

    // Looking through the reconstructed muon collection
    for (MAuint32 i=0;i<event.rec()->muons().size();i++)
    {
      const RecLeptonFormat& mu = event.rec()->muons()[i];
      cout << "----------------------------------" << endl;
      cout << "Muon" << endl;
      cout << "----------------------------------" << endl;
      cout << "index=" << i+1 
                << " charge=" << mu.charge() << endl;
      cout << "px=" << mu.px()
                << " py=" << mu.py()
                << " pz=" << mu.pz()
                << " e="  << mu.e()
                << " m="  << mu.m() << endl;
      cout << "pt=" << mu.pt() 
                << " eta=" << mu.eta() 
                << " phi=" << mu.phi() << endl;
      cout << "ET/PT isolation criterion =" << mu.ET_PT_isol() << endl;
      cout << "pointer address to the matching MC particle: " 
           << mu.mc() << endl;
      cout << endl;
    }

    // Looking through the reconstructed hadronic tau collection
    for (MAuint32 i=0;i<event.rec()->taus().size();i++)
    {
      const RecTauFormat& tau = event.rec()->taus()[i];
      cout << "----------------------------------" << endl;
      cout << "Tau" << endl;
      cout << "----------------------------------" << endl;
      cout << "tau: index=" << i+1 
                << " charge=" << tau.charge() << endl;
      cout << "px=" << tau.px()
                << " py=" << tau.py()
                << " pz=" << tau.pz()
                << " e="  << tau.e()
                << " m="  << tau.m() << endl;
      cout << "pt=" << tau.pt() 
                << " eta=" << tau.eta() 
                << " phi=" << tau.phi() << endl;
      cout << "pointer address to the matching MC particle: " 
           << tau.mc() << endl;
      cout << endl;
    }

    // Looking through the reconstructed jet collection
    for (MAuint32 i=0;i<event.rec()->jets().size();i++)
    {
      const RecJetFormat& jet = event.rec()->jets()[i];
      cout << "----------------------------------" << endl;
      cout << "Jet" << endl;
      cout << "----------------------------------" << endl;
      cout << "jet: index=" << i+1 
           << " charge=" << jet.charge() << endl;
      cout << "px=" << jet.px()
           << " py=" << jet.py()
           << " pz=" << jet.pz()
           << " e="  << jet.e()
           << " m="  << jet.m() << endl;
      cout << "pt=" << jet.pt() 
           << " eta=" << jet.eta() 
           << " phi=" << jet.phi() << endl;
      cout << "b-tag=" << jet.btag()
           << " true b-tag (before eventual efficiency)=" 
           << jet.true_btag() << endl;
      cout << "EE/HE=" << jet.EEoverHE()
           << " ntracks=" << jet.ntracks() << endl;
      cout << endl;
    }

    // Transverse missing energy (MET)
    cout << "MET pt=" << event.rec()->MET().pt()
         << " phi=" << event.rec()->MET().phi() << endl;
    cout << endl;

    // Transverse missing hadronic energy (MHT)
    cout << "MHT pt=" << event.rec()->MHT().pt()
              << " phi=" << event.rec()->MHT().phi() << endl;
    cout << endl;

    // Total transverse energy (TET)
    cout << "TET=" << event.rec()->TET() << endl;
    cout << endl;

    // Total transverse hadronic energy (THT)
    cout << "THT=" << event.rec()->THT() << endl;
    cout << endl;
  }
  */
  return true;
}

