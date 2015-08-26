
// EnergyPlus Headers
#include <EnergyPlus.hh>

// JSON Header
#include <json/json.h>

// GlobalSearch Header
#include <objective_function.hh>
#include <genetic_algorithm.hh>

namespace EnergyPlus {

namespace WindowASHRAE1588RP {

void
CreateASHRAE1588RPConstructions( int & ConstrNum, bool & ErrorsFound );

Json::Value read1588Database(std::string filePath);

void searchDatabaseKeysForInput(const std::string &constructionName, const std::string &fieldName, const std::vector< std::string > &keys, const std::string &input, bool &ErrorsFound);

void createDummyVariables();

void removeDummyVariables();

void calcWindowPerformance(Real64 T_in, Real64 T_out, Real64 v_ws, Real64 I_s);

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
	bool isContinuous(std::string trait);
	Real64 getTraitLowerBound(std::string trait);
	Real64 getTraitUpperBound(std::string trait);

private: // Data
	Json::Value db;

public: // Data
	Json::Value tests;
	std::vector <Real64> wavelengths;
	std::vector <Real64> solarSpectrum;
	std::vector <Real64> photopicResponse;

};

class FenSysObjFunc : public globalSearch::ObjectiveFunction
{
	public: // Creation
		FenSysObjFunc(
				const std::string& constructionName,
				const ASHRAE1588Database& database,
				const Real64& uFactorTarget,
				const Real64& shgcTarget);

	public: // Methods
		virtual Real64 call(const std::vector<Real64>& fenestrationTraits);

	private: // Data
		std::string constructionName_;
		ASHRAE1588Database database_;
		Real64 uFactorTarget_;
		Real64 shgcTarget_;
};

class FenestrationSystem {
public: // Creation
	FenestrationSystem(
		const std::string &constructionName,
		const ASHRAE1588Database &db,
		const Real64 &uFactorTarget,
		const Real64 &shgcTarget,
		std::vector <Real64> fenestrationTraits);

private: // Data
	ASHRAE1588Database database;
	std::string constructionName;

	// Trait indices
	int fenestrationTypeIndex;
	int numberOfPanesIndex;
	int glazingSubstrateIndex;
	int glazingCoatingIndex;
	int gasTypeIndex;
	int spacerTypeIndex;
	int frameTypeIndex;

	// Trait values or names
	std::string fenestrationType;
	int numberOfPanes;
	Real64 glazingThickness;
	std::string glazingSubstrateType;
	std::string glazingCoatingType;
	std::string gasType;
	Real64 gapThickness;
	std::string spacerType;
	std::string frameType;
	Real64 frameWidth;
	Real64 dividerWidth;

	// Output values
	Real64 fenestrationArea;
	Real64 fenestrationWidth;
	Real64 fenestrationHeight;
	Real64 glazingArea;
	Real64 frameArea;
	Real64 frameConductance;
	Real64 numHorizontalDividers;
	Real64 numVerticalDividers;

	Real64 uCOG;
	Real64 uEOG;
	Real64 frameEdgeRatio;


public: // Data
	Real64 uFactor;
	Real64 uFactorTarget;
	Real64 uFactorDiff;
	Real64 shgc;
	Real64 shgcTarget;
	Real64 shgcDiff;
	Real64 visibleTransmittance;
	Real64 error;
	bool matched;

public: // Methods
	void calculate();
	Json::Value generateOutput();


};

} // WindowASHRAE1588RP

} // EnergyPlus
