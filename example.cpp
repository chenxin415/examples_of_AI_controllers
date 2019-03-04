#include "./headers/propagate_intervals.h"
#include "../flowstar-template/Continuous.h"

using namespace std;
using namespace flowstar;

datatype offset_in_constraint_comb = constr_comb_offset;
// So the above data assumes, that this number is
// same for all the networks you have in the setting,
// and also all the networks are single output network
// with the first one giving the output

int main()
{
	/*
	 * Declaration of the state variables.
	 * The first one should always be the local time variable which can be
	 * viewed as a preserved variable. It is only used internally by the library.
	 */
	unsigned int numVars = 3;

	int x_id = stateVars.declareVar("x");
	int y_id = stateVars.declareVar("y");
	int u_id = stateVars.declareVar("u");

	int domainDim = numVars + 1;


	/*
	 * Define the first continuous dynamics.
	 */
	Expression_AST<Real> deriv_x("y - x^3", stateVars);
	Expression_AST<Real> deriv_y("u", stateVars);
	Expression_AST<Real> deriv_u("0", stateVars);

	vector<Expression_AST<Real> > ode_rhs(numVars);
	ode_rhs[x_id] = deriv_x;
	ode_rhs[y_id] = deriv_y;
	ode_rhs[u_id] = deriv_u;

	Deterministic_Continuous_Dynamics dynamics(ode_rhs);



	/*
	 * Specify the parameters for reachability computation.
	 */
	Computational_Setting setting;

	// stepsize and order
	setting.setFixedStepsize(0.02, 4);

	// time horizon
	setting.setTime(0.2);

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




	// Simple range propagation
	char controller_file[] = "../systems_with_networks/Ex_12/neural_network_controller" ;
	network_handler system_network(controller_file);



	/*
	 * Initial set can be a box which is represented by a vector of intervals.
	 * The i-th component denotes the initial set of the i-th state variable.
	 */
	Interval init_x(0.7,0.9), init_y(0.7,0.9), init_u;
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


	for(int k=0; k<30; ++k)
	{
		std::vector<Interval> NN_input;
		initial_set.intEvalNormal(NN_input, setting.tm_setting.step_exp_table, setting.tm_setting.cutoff_threshold);

		vector< vector< datatype > > input_interval(2, vector< datatype >(2,0));
		input_interval[0][0] = NN_input[0].inf();
		input_interval[0][1] = NN_input[0].sup();

		input_interval[1][0] = NN_input[1].inf();
		input_interval[1][1] = NN_input[1].sup();

		printf("[%lf, %lf], \t [%lf, %lf]\n", input_interval[0][0], input_interval[0][1],
				input_interval[1][0], input_interval[1][1]);

		vector< vector< datatype > > input_constraints;
		create_constraint_from_interval(input_constraints, input_interval);

		vector< datatype > output_range(2,0);
	 	system_network.return_interval_output(input_constraints, output_range, 1);

	 	cout << "output_range = [" << output_range[0] << " , " << output_range[1] << " ]" << endl;

	 	Interval u_range(output_range[0], output_range[1]);
	 	u_range -= 4;			// -4

	 	Real c;
	 	u_range.remove_midpoint(c);
	 	TaylorModel<Real> tm_u(c, domainDim);
	 	tm_u.remainder = u_range;

	 	initial_set.tmvPre.tms[2] = tm_u;

		printf("Step %d\n", k);


		int res = dynamics.reach(result, setting, initial_set, unsafeSet);

		if(res == COMPLETED_SAFE || res == COMPLETED_UNSAFE || res == COMPLETED_UNKNOWN)
		{
			initial_set = result.fp_end_of_time;
		}
		else
		{
			printf("Terminated due to too large overestimation.\n");
			break;
		}
	}


	// plot the flowpipes in the x-y plane
	result.transformToTaylorModels(setting);

	Plot_Setting plot_setting;
	plot_setting.setOutputDims(x_id, y_id);

	// you need to create a subdir named outputs
	plot_setting.plot_2D_interval_MATLAB("Ex_12_INT", result);

	return 0;
}
