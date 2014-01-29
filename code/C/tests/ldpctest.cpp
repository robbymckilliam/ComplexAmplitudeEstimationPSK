/* 
 * File:   ldpctest.cpp
 * Author: mckillrg
 *
 * Created on 29/01/2014, 1:58:37 PM
 */

#include <stdlib.h>
#include <iostream>

/*
 * Test Gottfried's LDPC code
 */

bool testencode() {
    
    
    
    bool pass = true;
    return pass;
}

void runtest(string name, function<bool()> test) {
    cout << name << " ... ";
    if (!test()) cout << "FAIL" << endl;
    else cout << "pass" << endl;
}

int main(int argc, char** argv) {
    runtest("test coherent Mackenthun estimate", testencode);
    return (EXIT_SUCCESS);
}

