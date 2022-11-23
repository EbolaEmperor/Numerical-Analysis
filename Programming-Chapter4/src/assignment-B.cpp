#include "FPN.h"
#include <iostream>

using std::cout;
using std::endl;

int main(){
    FPN_System fpns(-1, 1, 2, 3);
    FPN x(fpns);

    cout << "Question #1" << endl;
    x.turn_UFL();
    cout << "UFL: " << x.to_double() << endl;
    x.turn_OFL();
    cout << "OFL: " << x.to_double() << endl << endl;

    cout << "Question #2" << endl;
    x.turn_OFL(); x.turn_neg();
    int card = 0;
    do{
        card++;
        cout << x.to_double() << " ";
        if(card==30) exit(0);
    } while(x.turn_next());
    cout << endl << "Cardinality: " << card << endl << endl;

    cout << "Question #4" << endl;
    x.turn_UFL(); x.turn_neg();
    while(!x.isUFL()){
        card++;
        cout << x.to_double() << " ";
        x.turn_next_ext();
    }
    cout << endl << "Extended Cardinality: " << card << endl;

    return 0;
}