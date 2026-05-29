#include "SampleAnalyzer/User/Analyzer/emunu_pdata.h"
using namespace MA5;
using namespace std;
using namespace MA5;
using namespace std;
#include "TLorentzVector.h"

// Global storage for 3-momentum components
vector<double> px_muons, py_muons, pz_muons;
vector<double> px_electrons, py_electrons, pz_electrons;
vector<double> px_neutrinos, py_neutrinos;
vector<double> vis_mass;
//vector<double> coll_mass;
int nevent = 0;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool emunu_pdata::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  cout << "BEGIN Initialization" << endl;
  // initialize variables, histos
   // Clear vectors
    px_muons.clear(); py_muons.clear(); pz_muons.clear();
    px_electrons.clear(); py_electrons.clear(); pz_electrons.clear();
    px_neutrinos.clear(); py_neutrinos.clear(); vis_mass.clear();
  cout << "END   Initialization" << endl;
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void emunu_pdata::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
   cout << "BEGIN Finalization" << endl;
    ofstream myfile;
    myfile.open("output_momentum.txt");

    // Write momentum values to file
    size_t min_size = std::min(px_muons.size(), px_electrons.size());
	for (size_t i = 0; i < min_size; i++)
    {
        myfile << px_muons[i] << " " << py_muons[i] << " " << pz_muons[i] << " "
               << px_electrons[i] << " " << py_electrons[i] << " " << pz_electrons[i] << " "
               << px_neutrinos[i] << " " << py_neutrinos[i] << " " << vis_mass[i] << endl;
    }

    myfile.close();
    cout << "Saved " << px_electrons.size() << " events." << endl;
    cout << "END   Finalization" << endl;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool emunu_pdata::Execute(SampleFormat& sample, const EventFormat& event)
{
  double myWeight = 1.;
    if (!Configuration().IsNoEventWeight() && event.mc() != 0)
        myWeight = event.mc()->weight();

    Manager()->InitializeForNewEvent(myWeight);
	
	vector<TLorentzVector> jets;
	int njets = 0;
	int Nb = 0;


    if (event.rec() != 0)
    {
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

        // Extract electrons
        for (MAuint32 i = 0; i < event.rec()->electrons().size(); i++)
        {
            const RecLeptonFormat &elec = event.rec()->electrons()[i];
			if(elec.pt() > 5 && abs(elec.eta()) < 2.5)
			{
				ne5++;
			}
            if (elec.pt() > 10 && abs(elec.eta()) < 2.4)
            {
                ne++;
				TLorentzVector temp;
                temp.SetPxPyPzE(elec.px(), elec.py(), elec.pz(), elec.e());
                electrons.push_back(temp);
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
				njets++;
				TLorentzVector temp;
				temp.SetPxPyPzE(jet.px(),jet.py(),jet.pz(),jet.e());
				jets.push_back(temp);
			}
			if (jet.btag() >0)
			{
				Nb++;
			}
		}

        // Extract MET (neutrino 2-momentum)
		TLorentzVector met;
		met.SetXYZM	(event.rec()->MET().px(),event.rec()->MET().py(),0,0);
		
	

        // Store momenta only if we have exactly one electron and one muon, also have the correct delta phi
        if(ne == 1 && nmu == 1 && ne5 < 2 && nmu10 < 2 && ntau20 < 1 && njets == 1) //Modified for 1jet
		{
			if(electrons[0].DeltaR(muons[0]) > 0.3 && abs(electrons[0].DeltaPhi(met)) < 0.7 && abs(electrons[0].DeltaPhi(muons[0])) > 2.2)
			{
				TLorentzVector tavis = electrons[0];
				TLorentzVector mu_pt = muons[0];
				double mvis = (tavis+muons[0]).M();	
				double ptnuest = (tavis.Px()*met.Px()+tavis.Py()*met.Py())/tavis.Pt();
				double xtvis = tavis.Pt()/(tavis.Pt()+ptnuest);
				double computed_mcol = mvis/sqrt(xtvis);
				if (computed_mcol > 0)
					
				{
				// Store muon 3-momentum
				px_muons.push_back(muons[0].Px());
				py_muons.push_back(muons[0].Py());
				pz_muons.push_back(muons[0].Pz());

				// Store electron 3-momentum
				px_electrons.push_back(electrons[0].Px());
				py_electrons.push_back(electrons[0].Py());
				pz_electrons.push_back(electrons[0].Pz());

				// Store neutrino (MET) 2-momentum
				px_neutrinos.push_back(met.Px());
				py_neutrinos.push_back(met.Py());
				
				// Store Collinear Mass
				// Store Visible Mass
				vis_mass.push_back(mvis);	
				}
			}
		}
    }
    return true;
}

