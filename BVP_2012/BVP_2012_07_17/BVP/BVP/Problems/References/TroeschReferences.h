#pragma once
#include <map>
#include <initializer_list>
#include "../../Utils/AuxUtils.h"

/// <summary>
/// Reference data for the Troesch's problem
/// </summary>
class TroeschReferences
{
public:

	/// <summary>
	/// Data structure to represent reference data for the tangent
	/// of the solution in the lefr and right boundary points
	/// </summary>
	struct TangentReference
	{
		/// <summary>
		/// Constructor from the initializer list
		/// </summary>
		TangentReference(const std::initializer_list<double>& args) noexcept
		{
			if (args.size() == 0)
				return;

			int arg_id = 0;
			for (const auto val : args)
			{
				if (arg_id == 0)
					LeftPointTangent = val;

				if (arg_id == 1)
				{
					RightPointTangent = val;
					break;
				}

				arg_id++;
			}
		}

		/// <summary>
		/// Default constructor
		/// </summary>
		TangentReference() {}

		/// <summary>
		/// Left boundaty point rangent
		/// </summary>
		double LeftPointTangent = std::numeric_limits<double>::quiet_NaN();
		/// <summary>
		/// Right boundaty point tangent
		/// </summary>
		double RightPointTangent = std::numeric_limits<double>::quiet_NaN();
	};

	/// <summary>
	/// Reference data from
	/// @article{Gen_sol_of_TP,
	///author = "Vazquez-Leal, Hector and Khan, Yasir and Fern\'andez-Anaya, Guillermo
	///	and Herrera - May, Agust{ \'\i}n and Sarmiento-Reyes, Arturo and
	///	Filobello - Nino, Uriel and Jimenez - Fern\'andez, V{\'\i}ctor-M. and
	///	Pereyra - D{\'\i}az, Domitilo",
	///title = "{A general solution for Troesch's problem.}",
	///language = "English",
	///year = "2012",
	///doi = {10.1155 / 2012 / 208375},
	///journal = "Mathematical Problems in Engineering",
	///abstract = "{Summary: The homotopy perturbation method (HPM) is employed to obtain
	///	an approximate solution for the nonlinear differential equation which
	///	describes Troesch's problem. In contrast to other reported solutions
	///	obtained by using variational iteration method, decomposition method
	///	approximation, homotopy analysis method, Laplace transform decomposition
	///	method,and HPM method, the proposed solution shows the highest degree of
	///	accuracy in the results for a remarkable wide range of values of Troesch's
	///	parameter.}",
	///classmath = "{*65L99 (Numerical methods for ODE)
	///34A45(Theoretical approximation of solutions of ODE)
	///}"},
	/// </summary>
	const std::map<int, TangentReference> tangent_Vazquez_Leal = {
			{1, {8.45197542e-1, 1.341828780}},
			{2, {5.186322404e-1, 2.406790318}},
			{3, {2.55607567e-1, 4.266151411}},
			{4, {1.118803125e-1, 7.254574096}},
			{5, {4.575046433e-2, 1.210049478e+1}},
			{6, {1.795094954e-2, 2.003575791e+1}},
			{7, {6.867509691e-3, 3.308525498e+1}},
			{8, {2.587169418e-3, 5.457983465e+1}},
			{9, {9.655845408e-4, 9.000602248e+1}},
			{10, {3.583377845e-4, 1.484064126e+2}},
			{12, {4.891062176e-5, 4.034263503e+2}},
			{15, {2.444513025e-6, 1.808042673e+3}},
			{20, {1.648773182e-8, 2.202629966e+4}},
			{22, {2.231499933e-9}},
			{25, {1.111027228e-10}},
			{28, {5.531510890e-12}},
			{30, {7.486093793e-13}},
			{50, {1.542999878e-21}},
			{100, {2.976060781e-43}},
			{200, {1.107117221e-86}},
			{500, {5.699661125e-217}},
	};

	/// <summary>
	/// The following values of the solution slopes at the left and right boundary points
	/// have been calculated via Maple 2021 using this code
	/// Digits := 60;
	///l: = 'l';
	///error_base: = 0.1 * 10 ^ (-20);
	///err: = error_base;
	///for n to 30 do
	///	sol_control_c : = dsolve({ (1 - l) * (diff(u(x), x, x) - n ^ 2 * u(x)) + l * (diff(u(x), x, x) * (1 - u(x) ^ 2) + 2 * u(x) * diff(u(x), x) ^ 2 - n ^ 2 * u(x) * (1 + u(x) ^ 2)) = 0, u(0) = 0, u(1) = tanh(n / 4) }, numeric, output = listprocedure, maxmesh = 60000, continuation = l, range = 0 .. 1, abserr = err);
	/// sol_u_control_c: = rhs(sol_control_c[2]);
	///sol_u_prime_control_c: = rhs(sol_control_c[3]);
	///y: = x -> 4 * arctanh(sol_u_control_c(x)) / n;
	///y_prime: = x -> 4 * sol_u_prime_control_c(x) / (n * (1 - sol_u_control_c(x) ^ 2));
	///err: = error_base * y_prime(0);
	///print(n, evalf[20](y_prime(0)), evalf[20](y_prime(1)));
	///end do;
	/// </summary>
	const std::map<int, TangentReference> tangent_Maple_2021 = {
		{1, {0.84520268530995105992, 1.3418378623684903750}},
		{2, {0.51862121926934020856, 2.4069398312470712808}},
		{3, {0.25560421556293310992, 4.2662228618028235936}},
		{4, {0.11188016477074884192, 7.2545835747685823436}},
		{5, {0.045750461406318740100, 12.100495450777814326}},
		{6, {0.017950949489545842823, 20.035757896358684170}},
		{7, {0.0068675096950569237216, 33.085255288014834440}},
		{8, {0.0025871694189625792576, 54.579834455573441032}},
		{9, {0.00096558454107617375688, 90.006022309162966360}},
		{10, {0.00035833778463081369041, 148.40642115601013388}},
		{11, {0.00013252595622099004947, 244.68784549281781274}},
		{12, {0.000048910621759197887752, 403.42631474056142088}},
		{13, {0.000018028344597841978695, 665.14012960516910728}},
		{14, {0.0000066401087093429676108, 1096.6322465464930623}},
		{15, {0.0000024445130237432486686, 1808.0418613716930637}},
		{16, {8.9967757878636890664e-7, 2980.9576515791003722}},
		{17, {3.3106026949973454110e-7, 4914.7686368307653988}},
		{18, {1.2180976920741690444e-7, 8103.0838041655799428}},
		{19, {4.4815661929468488680e-8, 13359.726754810042546}},
		{20, {1.6487731827804035743e-8, 22026.465749406786042}},
		{21, {6.0657142766626446276e-9, 36315.502646710189174}},
		{22, {2.2314999333617785720e-9, 59874.141698496118468}},
		{23, {8.2093373809645267244e-10, 98715.771000630414592}},
		{24, {3.0200705232860191568e-10, 1.6275479141285965085e+5}},
		{25, {1.1110272283399722472e-10, 2.6833728651714794373e+5}},
		{26, {4.0872527453337473212e-11, 4.4241339200666022624e+5}},
		{27, {1.5036189304255441858e-11, 7.2941636984632986388e+5}},
		{28, {5.5315108863241889600e-12, 1.2026042841639472026e+6}},
		{29, {2.0349304652756557913e-12, 1.9827592635370682242e+6}},
		{30, {7.4860937950438118460e-13, 3.2690173724718024618e+6}},
		{31, {2.7539806648334038878e-13, 5.3896984762827488644e+6}},
		{32, {1.0131330159013064803e-13, 8.8861105205077436628e+6}},
		{33, {3.7271084072874085880e-14, 1.4650719428952882282e+7}},
		{34, {1.3711266317060284765e-14, 2.4154952753573807716e+7}},
		{35, {5.0440931548033120836e-15, 3.9824784397574635313e+7}},
		{36, {1.855618207672799e-15, 6.565996916237044e+7}},
		{37, {6.826437974477496E-16, 1.082549873554087E+8}},
		{38, {2.511306205497748E-16, 1.784823009553702E+8}},
		{39, {9.238579275622296e-17, 2.942675657296533e+8}},
		{40, {3.398683390222854E-17, 4.851651986914964E+8}},
		{41, {1.250305748341844E-17, 7.999021859604564E+8}},
		{42, {4.599617804459472e-18, 1.318815719945113e+9}},
		{43, {1.692104828516448E-18, 2.174359494205253E+9}},
		{44, {6.224905789434200E-19, 3.584912963483024E+9}},
		{45, {2.290014863664620e-19, 5.910521790684900e+9}},
		{46, {8.424493884314028E-20, 9.744803094508632E+9}},
		{47, {3.099198102563952E-20, 1.606646374522954E+10}},
		{48, {1.140131266106665e-20, 2.648911893017064e+10}},
	};

	/// <summary>
	/// Returns relative diference between the tro tangent values and the corresponding refetence values if they exists
	/// If the reference values does not exist the output is NaN
	/// </summary>
	TangentReference get_rel_diff(const int lambda, const double left_bpundary_tangent, const double right_boundary_tangent, const double maple_referencs = true)
	{
		const auto& data = maple_referencs ? tangent_Maple_2021 : tangent_Vazquez_Leal;

		const auto reference_ptr = data.find(lambda);

		if (reference_ptr == data.end())
			return TangentReference();

		const auto& reference = (*reference_ptr).second;

		return { auxutils::rel_diff(reference.LeftPointTangent, left_bpundary_tangent),
				auxutils::rel_diff(reference.RightPointTangent, right_boundary_tangent) };
	}
};