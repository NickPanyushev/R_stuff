unsigned long long rawNAN = 0x7ff8000000000000;
double NAN = *( double* ) &rawNAN;

unsigned long long rawPINF = 0x7ff0000000000000;
double PINF = *( double* ) &rawPINF;

unsigned long long rawNINF = 0xfff0000000000000;
double NINF = *( double* ) &rawNINF;
