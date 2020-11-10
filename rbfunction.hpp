#ifndef RBP
#define RBP
#include <eigen3/Eigen/Core>
/*This class takes care of the RBP functions of the meshless interpolator*/
class RBFunction{

	typedef double (*pRBF)(Eigen::VectorXd, Eigen::VectorXd, double c2);

    private:

        //stores the function if it is analytically defined
        pRBF function;
        
    public:

        //Initializes the class
        RBFunction(void);

        /*
        Initialization of the object with the function
        it is suposed to perform.
        input has to be of the kind:
        double (*input)(Eigen::VectorXd, Eigen::VectorXd, double);
        */
        void Init(pRBF input);

        /*
        Returns the value of the RBfunction in X with
        normal vector N.
        Inputs: Eigen::VectorXd, Eigen::VectorXd, double
        */
        double Value(Eigen::VectorXd x, 
                    Eigen::VectorXd x_i,
                    double c2);
};  

//Default function which always returns  0.0f
double Default_RBF(Eigen::VectorXd x, 
                  Eigen::VectorXd x_i,
                  double c2);
#endif