import os
import subprocess
import resource
import time

GENERATE_SOL = "GENERATE_SOL"
EXACT = "EXACT"

class Solver:
    def __init__(self, filename, args, mode):
        self.filename = filename
        self.args = args
        self.mode = mode
    def __str__(self):
        return "\"{} {}\" ({})".format(self.filename, " ".join(self.args), self.mode)


class Result:
    def __init__(self):
        self.error = None
        self.err_extra = ""
        self.opt = None
        self.expected_opt = 0

    def is_success(self):
        return self.error is None

    def __str__(self):
        out = ""
        if self.is_success():
            out = "SUCCESS"
        else:
            out = "ERROR({}): {}".format(self.error, self.err_extra)
        out = out + " Opt={} Expect={}".format(self.opt, self.expected_opt)
        return out
 
ERR_TIMEOUT = "timeout"
ERR_MEMORY = "memory"
ERR_RETCODE = "retcode"
ERR_PARSE = "parse"
ERR_WRONG = "wronganswer"
ERR_OTHER = "other"

def run_test(solver, test, timeout):
    print("Testing {} on {}".format(solver, test.filename))
    print("    ", end="")
    result = {
        "status": "ok",
        "file": test.filename
    }

    result["solver_cmdline"] = "_".join(solver.args)
    output_file = output_dir + "/" + test.filename[10:-2] + "_".join(solver.args) + ".result"

    if solver.mode == GENERATE_SOL:
        if os.path.exists(test.filename[:-2] + "ans"):
            print("answer file exists, not generating")
            return result
    else:
        if os.path.exists(output_file):
            result["status"] = "skip"
            return result
            
    start_cputime = resource.getrusage(resource.RUSAGE_CHILDREN)
    proc = subprocess.Popen([solver.filename, *solver.args, test.filename],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            #stderr=subprocess.PIPE,
            universal_newlines=True)
    try:
        stdout_data, stderr_data = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.terminate()
        stdout_data, stderr_data = proc.communicate()
        result["status"] = ERR_TIMEOUT
        result["wallclock_time"] = timeout * 1000
        end_cputime = resource.getrusage(resource.RUSAGE_CHILDREN)
        result["cpu_time"] = (end_cputime.ru_utime - start_cputime.ru_utime) * 1000
    proc.kill()

    if result["status"] != ERR_TIMEOUT and proc.returncode != 0:
        result["status"] = ERR_RETCODE
        result["error"] = "{}".format(proc.returncode)

    try:
        for line in stdout_data.splitlines():
            key, val = line.strip().split()
            result[key] = val
    except ValueError:
        result["status"] = ERR_PARSE
        result["error"] = stdout_data.strip()
        return result

    if result["status"] == "ok":
        try:
            with open(test.filename[:-4] + "ans", "r") as f:
                result["expected_sol_val"] = f.read().strip()
        except FileNotFoundError:
            print("answer file doesn't exist, generating")
            result["expected_sol_val"] = result["sol_val"].strip()
            with open(test.filename[:-4] + "ans", "w") as f:
                f.write(result["sol_val"])
                f.write("\n")

        if solver.mode == EXACT:
            if result["expected_sol_val"] != result["sol_val"]:
                result["status"] = ERR_WRONG
                result["error"] = "expected=={}".format(result["expected_sol_val"])

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w") as f:
        for key in result:
            f.write(str(key) + " " + str(result[key]) + "\n")
    
    return result

nfail = 0
def run_tests(solvers, tests, f, timeout=3):
    global nfail
    for solver in solvers:
        for test in tests:
            result = run_test(solver, test, timeout)
            print(result["status"], "cpu time: ", result["cpu_time"])
            if result["status"] not in ["ok", "skip"]:
                nfail += 1
            f.write(str(result) + "\n")
            f.flush()

class Test:
    def __init__(self, filename):
        self.filename = filename
        with open(filename, "r") as f:
            testdata = f.read()
        lines = testdata.split("\n")

    def __str__(self):
        return "{}".format(self.filename)

def get_all_tests():
    tests = []
    for group in sorted(os.listdir("instances")):
        files = os.listdir("instances/" + group)
        for file in sorted(files):
            if file[-5:] == ".msti":
                tests.append(Test("instances/" + group + "/" + file))
    return tests

output_dir = "test_results"
all_solvers = [
    Solver("build/mstisolver", ["-r-fraction", "0.05", "-dp-cost-scale", "100000", "-f"], EXACT),
    Solver("build/mstisolver", ["-r-fraction", "0.05", "-no-canonicalize", "-f"], EXACT),
    Solver("build/mstisolver", ["-q", "-f"], EXACT),
]

tests = get_all_tests()
wei2021_gen_tests = list(filter(lambda t: "wei2021integer_gen" in t.filename, tests))
wei2021_real_tests = list(filter(lambda t: "wei2021integer_real" in t.filename, tests))
small = list(filter(lambda t: "small" in t.filename, tests))
bazgan2012efficient = list(filter(lambda t: "bazgan2012efficient" in t.filename, tests))
new = list(filter(lambda t: "new" in t.filename, tests))
os.makedirs(output_dir, exist_ok=True)
with open(output_dir + "/all_results", "a") as f:
    run_tests([all_solvers[0]], wei2021_gen_tests, f, timeout=3600)
    run_tests([all_solvers[0]], wei2021_real_tests, f, timeout=3600)
    run_tests([all_solvers[2]], small, f, timeout=3600)
    # run_tests([all_solvers[2]], bazgan2012efficient, f, timeout=3600)
    # run_tests([all_solvers[2]], new, f, timeout=60)

print("\n{} FAILURES".format(nfail))
