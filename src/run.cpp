#include <bits/stdc++.h>
#include <filesystem>
#include <thread>
#include <mutex>
#include <chrono>

using namespace std;
namespace fs = std::filesystem;

// ===== main_2 선언 (QCSS 내부 함수) =====
int main_2(int argc, char** argv);

// =============================
// Config
// =============================
string dirpos   = "../example";
string binaryName = "./CCD";
string outpos   = "../output";
string logpos   = "../log";

mutex mtx;

// result: {benchname, [add_2q, fidelity, runtime]}
vector<pair<string, vector<string>>> benchresult;


// =========================================================
// Python run.py style: extract digits only (NO stoi)
// =========================================================
int extract_int(const string& s) {
    long long v = 0;
    bool found = false;
    for (char c : s) {
        if (isdigit(c)) {
            found = true;
            v = v * 10 + (c - '0');
        }
    }
    return found ? (int)v : 0;
}

string extract_numeric(const string& s) {
    string out;
    for (char c : s) {
        if (isdigit(c) || c == '.')
            out.push_back(c);
    }
    return out.empty() ? "0" : out;
}


// =========================================================
// execute external command + timeout (sh -c)
// =========================================================
bool runWithTimeoutAndTee(const string& cmd, const string& logPath, int time_limit_sec) {
    FILE* pipe = popen((cmd + " 2>&1").c_str(), "r");
    if (!pipe) return false;

    ofstream fout(logPath);
    if (!fout.is_open()) {
        pclose(pipe);
        return false;
    }

    auto start = chrono::steady_clock::now();

    char buffer[4096];
    while (true) {
        // read a line
        if (fgets(buffer, sizeof(buffer), pipe) == nullptr)
            break;

        // 1) print to terminal
        cout << buffer;

        // 2) write to logfile
        fout << buffer;

        // 3) check timeout
        auto now = chrono::steady_clock::now();
        if (chrono::duration_cast<chrono::seconds>(now - start).count() > time_limit_sec) {
            fout << "\n[TIMEOUT]\n";
            cout << "\n[TIMEOUT]\n";
            fout.close();
            pclose(pipe);
            return false;
        }
    }

    fout.close();
    pclose(pipe);
    return true;
}

// =========================================================
// LOG PARSER
// =========================================================
vector<string> parseLog(const string& path) {
    ifstream f(path);
    if (!f.is_open()) return {"0","0","0"};

    string line;
    int add_2q = 0;
    int fidelity = 0;
    string runtime = "0";

    while (getline(f, line)) {
        if (line.find("add_2q") != string::npos) {
            add_2q = extract_int(line.substr(line.find(":")+1));
        }
        if (line.find("fidelity") != string::npos) {
            fidelity = extract_int(line.substr(line.find(":")+1));
        }
        if (line.find("END PROGRAM") != string::npos) {
            string tok, last;
            stringstream ss(line);
            while (ss >> tok) last = tok;
            runtime = extract_numeric(last);
        }
    }
    return {to_string(add_2q), to_string(fidelity), runtime};
}


// =========================================================
// execute one benchmark
// =========================================================
void executeBench(const string& bench, int timeout) {
    string input  = dirpos + "/" + bench;
    string output = outpos + "/" + bench;

    char buf[32];
    time_t t = time(nullptr);
    strftime(buf, 32, "%m_%d_%H_%M", localtime(&t));
    string curtime = buf;

    string logfile = logpos + "/" + bench + "." + curtime + ".log";

    //  tee 없이 순수한 명령만 실행 (tee는 C++에서 직접 구현)
    string exe = binaryName + " " + input + " " + output;

    cout << "\n=== Running: " << exe << " ===\n";

    bool ok = runWithTimeoutAndTee(exe, logfile, timeout);

    if (!ok) {
        cout << "[TIMEOUT] " << bench << "\n";
        return;
    }

    vector<string> parsed = parseLog(logfile);

    {
        lock_guard<mutex> lk(mtx);
        string name = bench.substr(0, bench.size() - 5);
        benchresult.push_back({name, parsed});
    }

    cout << "Successfully run!\n";
}



// =========================================================
// MAIN
// =========================================================
int main(int argc, char* argv[]) {

    // ---------------------------
    // 0) 인자 검사
    // ---------------------------
    if (argc <= 1) {
        cout << "usage:\n"
             << "  QCSS MODE: ./Convergent_Capstone_Design input.qasm output.qasm\n"
             << "  RUN MODE : ./Convergent_Capstone_Design all | list | index\n";
        return 0;
    }

    string arg1 = argv[1];

    // ---------------------------
    // 1) QCSS MODE 판별
    // QCSS 모드: arg1이 ".qasm"으로 끝나면 그냥 main_2 실행
    // ---------------------------
    if (arg1.size() > 5 && arg1.substr(arg1.size()-5) == ".qasm") {

        // QCSS pipeline 실행
        return main_2(argc, argv);
    }

    // ---------------------------
    // 2) RUN MODE
    // run.py 포팅 기능 수행
    // ---------------------------

    vector<string> benchlist;

    for (auto& p : fs::directory_iterator(dirpos)) {
        if (p.path().extension() == ".qasm") {
            benchlist.push_back(p.path().filename().string());
        }
    }
    sort(benchlist.begin(), benchlist.end());

    vector<string> toRun;

    // list
    if (arg1 == "list") {
        for (int i = 0; i < benchlist.size(); i++)
            cout << i << " : " << benchlist[i] << "\n";
        return 0;
    }

    // all
    if (arg1 == "all") {
        toRun = benchlist;
    }
    // index mode
    else if (all_of(arg1.begin(), arg1.end(), ::isdigit)) {
        int idx = stoi(arg1);
        if (idx >= 0 && idx < benchlist.size())
            toRun.push_back(benchlist[idx]);
    }
    // multiple index
    else if (argc >= 3 && all_of(argv[1], argv[argc-1], ::isdigit)) {
        for (int i = 1; i < argc; i++) {
            int idx = stoi(argv[i]);
            if (idx >= 0 && idx < benchlist.size())
                toRun.push_back(benchlist[idx]);
        }
    }
    else {
        // bench name
        toRun.push_back(arg1);
    }

    // ----- run benches -----
    vector<thread> threads;
    for (auto& b : toRun)
        threads.emplace_back(executeBench, b, 10000);

    for (auto& th : threads) th.join();

    sort(benchresult.begin(), benchresult.end(),
         [](auto& a, auto& b){ return a.first < b.first; });

    // ---- print result ----
    cout << "\n\n";
    cout << "----------------------------- Result -------------------------------\n";
    cout << "                        |                   Our                     |\n";
    cout << "      bench_name        | add_2q_num  |    fidelity    |   RunTime  |\n";
    cout << "--------------------------------------------------------------------\n";

    for (auto& br : benchresult) {
        if(br.second[0].find("2147483647") != string::npos) // timeout
        {cout << left << setw(23) << br.first << " | "
                 << right << setw(11) << "OOB" << " | "
                 << setw(14) << br.second[1] << " | "
                 << setw(10) << br.second[2] << " |\n";}
            
        else
           {cout << left << setw(23) << br.first << " | "
                 << right << setw(11) << br.second[0] << " | "
                 << setw(14) << br.second[1] << " | "
                 << setw(10) << br.second[2] << " |\n";} 
        
    }
    cout << "--------------------------------------------------------------------\n";

    return 0;
}
