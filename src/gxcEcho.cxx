#include "gxcEcho.h"

void EchoStr(const std::string &Str2Echo)
{
    std::cout << Str2Echo << std::endl;
    std::fstream logFile("gxmc.log", std::fstream::app);
    logFile << Str2Echo << std::endl;
    logFile.close();
}

void EchoMat(const std::string VarName, const std::vector<std::vector<gxcData>> &Var2Echo)
{
    int m = Var2Echo.size();
    int n = Var2Echo[0].size();

    std::cout << VarName << "=" << std::endl;

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
            std::cout << Var2Echo[i][j] << " ";
        std::cout << std::endl;
    }
}

void LogInitiate()
{
    std::fstream logFile("gxmc.log", std::fstream::in | std::fstream::out | std::ofstream::trunc);
    auto startTime = std::chrono::system_clock::now();
    std::time_t outputTime = std::chrono::system_clock::to_time_t(startTime);

    std::cout << "The simulation starts at " << std::ctime(&outputTime) << std::endl;

    logFile << "******** GXMC log file ********" << std::endl
            << "The simulation starts at " << std::ctime(&outputTime) << std::endl;

    logFile.close();
}
void EchoError(const std::string &Str2Echo)
{
    std::cout << "ERROR!!! "<<Str2Echo << std::endl;
    std::fstream logFile("gxmc.log", std::fstream::app);
    logFile << "ERROR!!! "<<Str2Echo << std::endl;
    logFile.close();
    exit(1);
}