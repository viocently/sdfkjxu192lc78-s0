#ifndef __LOG_H__
#define __LOG_H__
#include<chrono>
#include<iostream>
#include<fstream>
#include<string>
#include<malloc.h>

using namespace std;

typedef bool TIMETYPE;

const int DATEFLAG = false;
const int TIMEFLAG = true;

inline string getCurrentSystemTime(TIMETYPE type)
{
    auto tt = chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    struct tm* ptm = localtime(&tt);
    char date[60] = { 0 };
    if ( type == DATEFLAG )
        sprintf( date, "%d-%02d-%02d", (int)ptm->tm_year + 1900, (int)ptm->tm_mon + 1, (int)ptm->tm_mday );
    else 
        sprintf( date, "%02d:%02d:%02d", (int) ptm->tm_hour, (int)ptm->tm_min, (int)ptm -> tm_sec);

    return string(date);
}

inline void logger( string logMsg )
{
    string filePath = "LOG/log_" + getCurrentSystemTime( DATEFLAG ) + ".log"; 
    string now = getCurrentSystemTime(TIMEFLAG);
    ofstream os ( filePath.c_str(), ios_base::out | ios_base::app ); 
    os << now << "\t" << logMsg << "\n";
    os.close();
}



#ifndef _WIN32
#include<unistd.h>
inline void showProcessMemUsage()
{
    double virMem = 0;
    double resMem = 0;

    fstream fs;
    fs.open("/proc/self/stat", ios_base::in);

    int dataCnt = 0;
    string dataStr;
    unsigned long virMemB, resMemP;
    while (fs >> dataStr)
    {
        dataCnt++;
        if (dataCnt == 22)
            break;
    }
    fs >> virMemB >> resMemP;

    fs.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
    virMem = virMemB / 1024.0 / 1024.0 / 1024.0;
    resMem = resMemP * page_size_kb / 1024.0 / 1024.0;

    logger("Current virtual memory usage : " + to_string(virMem) + "G.");
    logger("Current resident memory usage : " + to_string(resMem) + "G.");
}
#endif

#endif
