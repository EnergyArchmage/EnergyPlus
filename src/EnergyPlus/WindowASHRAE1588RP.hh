
// EnergyPlus Headers
#include <EnergyPlus.hh>

// JSON Header
#include <json/json.h>

namespace EnergyPlus {

namespace WindowASHRAE1588RP {

void
CreateASHRAE1588RPConstructions( int & ConstrNum, bool & ErrorsFound );

Json::Value read_1588_database(std::string file_path);

void search_database_keys_for_input(const std::string &construction_name, const std::string &field_name, const std::vector< std::string > &keys, const std::string &input, bool &ErrorsFound);

void create_dummy_variables();

void remove_dummy_variables();

void calc_window_performance(Real64 T_in, Real64 T_out, Real64 v_ws, Real64 I_s);

class ASHRAE1588Database {
public: // Creation
	ASHRAE1588Database(Json::Value db);

public: // Methods
	std::vector <std::string> getTraitOptions(std::string trait);
	int getTraitIndexWithMaxUtility(std::string trait);
	std::string getTraitNameWithMaxUtility(std::string trait);
	Json::Value getDiscreteTraitValueWithMaxUtility(std::string trait);
	Real64 getContinuousTraitValueWithMaxUtility(std::string trait);
	int getTraitIndexByName(std::string trait, std::string name);
	Json::Value getTraitValueByName(std::string trait, std::string name);
	std::string getTraitNameByIndex(std::string trait, int index);
	Json::Value getTraitValueByIndex(std::string trait, int index);

private: // Data
	Json::Value db;

public: // Data
	Json::Value tests;
	std::vector <Real64> wavelengths;
	std::vector <Real64> solarSpectrum;
	std::vector <Real64> photopicResponse;

};

} // WindowASHRAE1588RP

} // EnergyPlus
