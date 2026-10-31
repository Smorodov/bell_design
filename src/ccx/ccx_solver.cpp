#include "ccx_solver.hpp"
#include "ccx_input_writer.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <utility>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <cerrno>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace BellOptimizer {
namespace {

bool valid_job_name(const std::string& job) {
    if (job.empty() || job == "." || job == "..") return false;
    return std::all_of(job.begin(), job.end(), [](unsigned char c) {
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
               (c >= '0' && c <= '9') || c == '_' || c == '-';
    });
}

#ifdef _WIN32
std::string quote_windows_argument(const std::string& arg) {
    std::string out = "\"";
    std::size_t slashes = 0;
    for (char c : arg) {
        if (c == '\\') {
            ++slashes;
        } else if (c == '"') {
            out.append(slashes * 2 + 1, '\\');
            out += '"';
            slashes = 0;
        } else {
            out.append(slashes, '\\');
            slashes = 0;
            out += c;
        }
    }
    out.append(slashes * 2, '\\');
    out += '"';
    return out;
}

bool run_process(const std::string& executable, const std::string& job,
                 const std::filesystem::path& work, bool quiet) {
    STARTUPINFOA startup{};
    startup.cb = sizeof(startup);
    PROCESS_INFORMATION process{};
    std::string command_line = quote_windows_argument(executable) + " " + quote_windows_argument(job);
    HANDLE null_handle = INVALID_HANDLE_VALUE;
    SECURITY_ATTRIBUTES sa{};
    sa.nLength = sizeof(sa);
    sa.bInheritHandle = TRUE;
    if (quiet) {
        null_handle = CreateFileA("NUL", GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE,
                                  &sa, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
        if (null_handle == INVALID_HANDLE_VALUE) return false;
        startup.dwFlags |= STARTF_USESTDHANDLES;
        startup.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
        startup.hStdOutput = null_handle;
        startup.hStdError = null_handle;
    }
    const BOOL created = CreateProcessA(executable.c_str(), command_line.data(), nullptr, nullptr,
                                        quiet ? TRUE : FALSE, 0, nullptr, work.string().c_str(),
                                        &startup, &process);
    if (null_handle != INVALID_HANDLE_VALUE) CloseHandle(null_handle);
    if (!created) return false;
    WaitForSingleObject(process.hProcess, INFINITE);
    DWORD exit_code = 1;
    const BOOL got_code = GetExitCodeProcess(process.hProcess, &exit_code);
    CloseHandle(process.hThread);
    CloseHandle(process.hProcess);
    return got_code && exit_code == 0;
}
#else
bool run_process(const std::string& executable, const std::string& job,
                 const std::filesystem::path& work, bool quiet) {
    const pid_t pid = fork();
    if (pid < 0) return false;
    if (pid == 0) {
        if (chdir(work.c_str()) != 0) _exit(126);
        if (quiet) {
            const int null_fd = open("/dev/null", O_WRONLY);
            if (null_fd < 0) _exit(126);
            if (dup2(null_fd, STDOUT_FILENO) < 0 || dup2(null_fd, STDERR_FILENO) < 0) _exit(126);
            if (null_fd > STDERR_FILENO) close(null_fd);
        }
        if (executable.find('/') == std::string::npos)
            execlp(executable.c_str(), executable.c_str(), job.c_str(), static_cast<char*>(nullptr));
        else
            execl(executable.c_str(), executable.c_str(), job.c_str(), static_cast<char*>(nullptr));
        _exit(127);
    }
    int status = 0;
    while (waitpid(pid, &status, 0) < 0) {
        if (errno == EINTR) continue;
        return false;
    }
    return WIFEXITED(status) && WEXITSTATUS(status) == 0;
}
#endif

} // namespace

CalculiXSolver::CalculiXSolver(std::string executable, bool quiet)
    : executable_(std::move(executable)), quiet_(quiet) {}

bool CalculiXSolver::solve(const ModalProblem& problem, ModalData& result) {
    result.clear();
    if (executable_.empty()) {
        if (!quiet_) std::cerr << "CalculiX executable path is empty\n";
        return false;
    }

    const std::filesystem::path work = std::filesystem::absolute(problem.work_directory);
    if (!std::filesystem::is_directory(work)) {
        if (!quiet_) std::cerr << "CalculiX work directory not found: " << work << '\n';
        return false;
    }
    const std::string job = problem.job_name.empty() ? "bell" : problem.job_name;
    if (!valid_job_name(job)) {
        if (!quiet_) std::cerr << "CalculiX: invalid job name (allowed: letters, digits, '_' and '-')\n";
        return false;
    }
    const std::filesystem::path input = work / (job + ".inp");
    const std::filesystem::path output = work / (job + ".dat");

    if (!CCXInputWriter::write(problem.mesh, problem.material, input.string())) {
        if (!quiet_) std::cerr << "CalculiX: failed to write input " << input << '\n';
        return false;
    }

    if (!run_process(executable_, job, work, quiet_)) {
        if (!quiet_) std::cerr << "CalculiX failed to launch or returned a non-zero exit code\n";
        return false;
    }

    result = CCXParser::parse_file(output.string(), quiet_);
    result.source_file = output.string();
    return !result.modes.empty();
}

} // namespace BellOptimizer
