#include "../flowstar/Continuous.h"
#include "../Bernstein_Polynomial_Approximation/bernstein_poly_approx.h"

using namespace std;
using namespace flowstar;


int main()
{
	// Declaration of the state variables.
	unsigned int numVars = 3;

	int d_err_id = stateVars.declareVar("d_err");
	int theta_err_id = stateVars.declareVar("t_err");
	int u_id = stateVars.declareVar("u");

	int domainDim = numVars + 1;


	// Define the continuous dynamics.
//	Expression_AST<Real> deriv_d_err("-sin(0.62831853071 - t_err)*0.58778525229 + cos(0.62831853071 - t_err)*0.80901699437");   // theta_r = pi/5
	Expression_AST<Real> deriv_d_err("cos(0 - t_err)");  // theta_r = 0
	Expression_AST<Real> deriv_theta_err("-u");
	Expression_AST<Real> deriv_u("0");

	vector<Expression_AST<Real> > ode_rhs(numVars);
	ode_rhs[d_err_id] = deriv_d_err;
	ode_rhs[theta_err_id] = deriv_theta_err;
	ode_rhs[u_id] = deriv_u;

	Deterministic_Continuous_Dynamics dynamics(ode_rhs);



	
	// Specify the parameters for reachability computation.
	
	Computational_Setting setting;

	unsigned int order = 4;

	// stepsize and order for reachability analysis
	setting.setFixedStepsize(0.02, order);

	// time horizon for a single control step
	setting.setTime(1);

	// cutoff threshold
	setting.setCutoffThreshold(1e-10);

	// queue size for the symbolic remainder
	setting.setQueueSize(200);

	// print out the steps
	setting.printOn();

	// remainder estimation
	Interval I(-0.01, 0.01);
	vector<Interval> remainder_estimation(numVars, I);
	setting.setRemainderEstimation(remainder_estimation);

	setting.prepare();

	/*
	 * Initial set can be a box which is represented by a vector of intervals.
	 * The i-th component denotes the initial set of the i-th state variable.
	 */
	Interval init_x(0,0.01), init_y(0,0.01), init_u(0);
	std::vector<Interval> X0;
	X0.push_back(init_x);
	X0.push_back(init_y);
	X0.push_back(init_u);




	// translate the initial set to a flowpipe
	Flowpipe initial_set(X0);

	// no unsafe set
	vector<Constraint> unsafeSet;

	// result of the reachability computation
	Result_of_Reachability result;

	// define the neural network controller
	char const *module_name = "dubins_controller_poly_approx";
	char const *function_name1 = "dubins_poly_controller";
	char const *function_name2 = "poly_approx_error";
	char const *degree_bound = "[3, 3]";
	char const *lips = "2.011";

	// perform 20 control steps
	for(int iter=0; iter<20; ++iter)
	{
		vector<Interval> box;
		initial_set.intEval(box, order, setting.tm_setting.cutoff_threshold);

		string strBox = "[" + box[0].toString() + "," + box[1].toString() + "]";

	
		string strExpU = bernsteinPolyApproximation(module_name, function_name1, degree_bound, strBox.c_str(), lips);
		double err = stod(bernsteinPolyApproximation(module_name, function_name2, degree_bound, strBox.c_str(), lips));

		Expression_AST<Real> exp_u(strExpU);

		TaylorModel<Real> tm_u;
		exp_u.evaluate(tm_u, initial_set.tmvPre.tms, order, initial_set.domain, setting.tm_setting.cutoff_threshold, setting.g_setting);

		tm_u.remainder.bloat(err);
	
		initial_set.tmvPre.tms[u_id] = tm_u;

		dynamics.reach(result, setting, initial_set, unsafeSet);

		if(result.status == COMPLETED_SAFE || result.status == COMPLETED_UNSAFE || result.status == COMPLETED_UNKNOWN)
		{
			initial_set = result.fp_end_of_time;
		}
		else
		{
			printf("Terminated due to too large overestimation.\n");
		}
	}


	// plot the flowpipes in the x-y plane
	result.transformToTaylorModels(setting);

	Plot_Setting plot_setting;
	plot_setting.setOutputDims(d_err_id, theta_err_id);

	mkres = mkdir("./outputs", S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for images.\n");
		exit(1);
	}

	// you need to create a subdir named outputs
	// the file name is example.m and it is put in the subdir outputs
	plot_setting.plot_2D_interval_MATLAB("example", result);

	return 0;
}
