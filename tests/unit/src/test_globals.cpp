/**
 * @file test_globals.cpp
 * @brief Unit tests for globals.cpp: parse_args, set_thread_ram_steps, parent_path,
 *        create_directory, init_timers, the print_* formatters, and the string tables.
 */
#include <unistd.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "../../../src/defs.h"
#include "../../../src/globals.h"
#include "../../../src/timer.h"
#include "test_helpers.h"

namespace {

/** @brief Returns the names of the given timers, in order. */
std::vector<std::string> timer_names(std::vector<timer> & timers) {
    std::vector<std::string> names;
    for (timer & t : timers) names.push_back(t.get_name());
    return names;
}

/* parent_path ************************************************************************************/

TEST(ParentPath, Nested) {
    EXPECT_EQ("a/b/", parent_path("a/b/c"));
}

TEST(ParentPath, SingleDir) {
    EXPECT_EQ("out/", parent_path("out/prefix"));
}

TEST(ParentPath, BareFilename) {
    // no separator anywhere, so there is no parent to return
    EXPECT_EQ("", parent_path("prefix"));
}

TEST(ParentPath, TrailingSlash) {
    // the trailing separator is itself the last one found, so the path is returned unchanged
    EXPECT_EQ("a/b/", parent_path("a/b/"));
}

TEST(ParentPath, Root) {
    EXPECT_EQ("/", parent_path("/"));
}

TEST(ParentPath, AbsoluteFile) {
    EXPECT_EQ("/", parent_path("/prefix"));
}

TEST(ParentPath, Empty) {
    EXPECT_EQ("", parent_path(""));
}

TEST(ParentPath, DotSlash) {
    EXPECT_EQ("./", parent_path("./prefix"));
}

TEST(ParentPath, DotDot) {
    EXPECT_EQ("../out/", parent_path("../out/prefix"));
}

/* create_directory *******************************************************************************/

TEST(CreateDirectory, SingleDir) {
    TempDir dir;
    create_directory(dir.path("a/b"));

    // only components followed by a separator are created, so the final "b" is not
    EXPECT_TRUE(std::filesystem::is_directory(dir.path("a")));
    EXPECT_FALSE(std::filesystem::exists(dir.path("a/b")));
}

TEST(CreateDirectory, NestedTrailingSlash) {
    TempDir dir;
    create_directory(dir.path("x/y/"));

    // the trailing separator makes "y" an intermediate component, so every level is created
    EXPECT_TRUE(std::filesystem::is_directory(dir.path("x")));
    EXPECT_TRUE(std::filesystem::is_directory(dir.path("x/y")));
}

TEST(CreateDirectory, AlreadyExistsOk) {
    TempDir dir;

    // every component of the first call already exists on the second, and EEXIST is tolerated
    create_directory(dir.path("dup/"));
    create_directory(dir.path("dup/"));
    EXPECT_TRUE(std::filesystem::is_directory(dir.path("dup")));
}

TEST(CreateDirectory, AbsoluteSkipsRoot) {
    // the scan starts one character past the front of the path, so a leading separator is never
    // treated as one: a single top-level component yields no mkdir call at all
    create_directory("/vcfdist_create_directory_should_be_a_noop");
    EXPECT_FALSE(std::filesystem::exists("/vcfdist_create_directory_should_be_a_noop"));

    // and a deeper absolute path succeeds without ever attempting to create the root
    TempDir dir;
    create_directory(dir.path("abs/leaf"));
    EXPECT_TRUE(std::filesystem::is_directory(dir.path("abs")));
}

TEST(CreateDirectory, MkdirFailureErrors) {
    GlobalsGuard guard;
    TempDir dir;

    // a regular file where a directory is needed makes mkdir fail with ENOTDIR, not EEXIST
    std::ofstream(dir.path("file")) << "not a directory\n";
    ASSERT_TRUE(std::filesystem::is_regular_file(dir.path("file")));

    EXPECT_EXIT(create_directory(dir.path("file/sub/")), testing::ExitedWithCode(1),
            "Unable to create directory");
}

TEST(CreateDirectory, NoSlashNoop) {
    const std::string relative = "vcfdist_create_directory_no_slash_noop";
    ASSERT_FALSE(std::filesystem::exists(relative));

    // with no separator there is nothing to iterate over, so no directory is created
    create_directory(relative);
    EXPECT_FALSE(std::filesystem::exists(relative));
}

TEST(CreateDirectory, EmptyNoop) {
    TempDir dir;
    const std::filesystem::path before = std::filesystem::current_path();
    std::filesystem::current_path(dir.path());

    // an empty path names no directory, so it returns before touching the filesystem
    create_directory("");
    EXPECT_TRUE(std::filesystem::is_empty(dir.path()));

    std::filesystem::current_path(before);
}

/* init_timers ************************************************************************************/

TEST(InitTimers, Populates) {
    GlobalsGuard guard;
    g.timers.clear();

    g.init_timers(timer_strs);

    ASSERT_EQ(idx(TIME_TOTAL)+1, g.timers.size());
    EXPECT_EQ(timer_strs, timer_names(g.timers));
    EXPECT_EQ("reading", g.stage(TIME_READ).get_name());
    EXPECT_EQ("total", g.stage(TIME_TOTAL).get_name());
}

TEST(InitTimers, EmptyInput) {
    GlobalsGuard guard;
    g.timers.clear();

    g.init_timers({});

    EXPECT_TRUE(g.timers.empty());
}

TEST(InitTimers, AppendsNotClears) {
    GlobalsGuard guard;
    g.timers.clear();

    // each call pushes onto the existing vector rather than replacing it
    g.init_timers({"first"});
    g.init_timers({"second", "third"});

    ASSERT_EQ(size_t(3), g.timers.size());
    EXPECT_EQ(std::vector<std::string>({"first", "second", "third"}), timer_names(g.timers));
}

TEST(InitTimers, WritesThisNotGlobal) {
    GlobalsGuard guard;
    g.timers.clear();

    // the timers land on the instance the method was called on, not on the global `g`
    Globals local;
    local.init_timers({"local"});

    ASSERT_EQ(size_t(1), local.timers.size());
    EXPECT_EQ("local", local.timers[0].get_name());
    EXPECT_TRUE(g.timers.empty());
}

/* String lookup tables ***************************************************************************/

TEST(StringTables, SizesMatchCount) {
    // these five tables are sized exactly to their governing count constant
    EXPECT_EQ(size_t(TYPES), type_strs.size());
    EXPECT_EQ(size_t(VARTYPES), vartype_strs.size());
    EXPECT_EQ(size_t(ERRTYPES), error_strs.size());
    EXPECT_EQ(size_t(SWITCHTYPES), switch_strs.size());
    EXPECT_EQ(size_t(CALLSETS), callset_strs.size());
}

TEST(StringTables, SizesWithSentinel) {
    // these four are one longer than their count constant, because a sentinel value is also a
    // valid subscript; shortening any of them to its count would read out of bounds
    EXPECT_EQ(AC_ERRTYPE_SLOTS, ac_strs.size());
    EXPECT_EQ(size_t(PHASES+1), phase_strs.size());
    EXPECT_EQ(".", ac_strs[AC_UNKNOWN]);
    EXPECT_EQ(".", phase_strs[PHASE_NONE]);

    // gt_strs and timer_strs have no count constant, so the highest valid index bounds them
    EXPECT_EQ(size_t(GT_OTHER+1), gt_strs.size());
    EXPECT_EQ(idx(TIME_TOTAL)+1, timer_strs.size());
    EXPECT_EQ("X|Y", gt_strs[GT_OTHER]);
    EXPECT_EQ("total", timer_strs[idx(TIME_TOTAL)]);
}

TEST(StringTables, IndexMapping) {
    EXPECT_EQ("TP", error_strs[ERRTYPE_TP]);
    EXPECT_EQ("1|1", gt_strs[GT_ALT1_ALT1]);
    EXPECT_EQ(".|.", gt_strs[GT_MISSING]);
    EXPECT_EQ("X|.", gt_strs[GT_HALF]);
    EXPECT_EQ("SNP", type_strs[TYPE_SUB]);
    EXPECT_EQ("QUERY", callset_strs[QUERY]);
    EXPECT_EQ("TRUTH", callset_strs[TRUTH]);
    EXPECT_EQ("SV", vartype_strs[VARTYPE_SV]);
    EXPECT_EQ("SWITCH", switch_strs[SWITCHTYPE_SWITCH]);
    EXPECT_EQ("INSIDE", region_strs[BED_INSIDE]);
    EXPECT_EQ("REF", type_strs[TYPE_REF]);
    EXPECT_EQ("CPX", type_strs[TYPE_CPX]);
    EXPECT_EQ("??", error_strs[ERRTYPE_UN]);
}

TEST(StringTables, AcSparse) {
    // only the two transitions that change the allele count while still counting as a TP carry a
    // non-"." string; every other allele-count error type prints as "."
    EXPECT_EQ("+", ac_strs[AC_ERR_1_TO_2]);
    EXPECT_EQ("-", ac_strs[AC_ERR_2_TO_1]);
    for (ac_errtype_t ac : EnumRange<ac_errtype_t, AC_ERRTYPE_SLOTS>{}) {
        if (ac == AC_ERR_1_TO_2 || ac == AC_ERR_2_TO_1) continue;
        EXPECT_EQ(".", ac_strs[ac]) << "ac_strs[" << idx(ac) << "]";
    }
}

TEST(StringTables, RegionUnpadded) {
    // these reach the LOCATION column of query.tsv and truth.tsv verbatim, so they must carry no
    // whitespace: an exact string comparison against the documented domain has to match
    EXPECT_EQ("OUTSIDE", region_strs[BED_OUTSIDE]);
    EXPECT_EQ("INSIDE", region_strs[BED_INSIDE]);
    EXPECT_EQ("BORDER", region_strs[BED_BORDER]);
    EXPECT_EQ("OFF_CTG", region_strs[BED_OFFCTG]);

    ASSERT_EQ(BEDLOC_SLOTS, region_strs.size());
    for (const std::string & s : region_strs) {
        EXPECT_EQ(std::string::npos, s.find(' ')) << "padded entry '" << s << "'";
    }
}

/* parse_args *************************************************************************************/

/**
 * @class ArgsFixture
 * @brief Openable files for the three mandatory arguments, plus a full-argv builder.
 *
 * parse_args() opens argv[1] and argv[2] with bcf_open() and argv[3] with fopen(), so every case
 * that reaches the optional-argument loop needs three files that exist. The temporary directory
 * owns the two written VCFs and doubles as scratch space for cases that write output.
 */
class ArgsFixture {
public:
    /** @brief Writes a header-only query and truth VCF, and points the reference at tiny.fasta. */
    ArgsFixture() {
        vcf_opts query_opts = make_vcf_opts(QUERY);
        query_opts.filename = "query.vcf";
        this->query_fn = write_tmp_vcf(this->dir, {}, query_opts);
        vcf_opts truth_opts = make_vcf_opts(TRUTH);
        truth_opts.filename = "truth.vcf";
        this->truth_fn = write_tmp_vcf(this->dir, {}, truth_opts);
        // data_path() resolves relative to the working directory, which the relative-prefix cases
        // change out from under the parse
        this->ref_fn = std::filesystem::absolute(data_path("tiny.fasta")).string();
    }

    /** @brief Returns the program name, the three mandatory arguments, then the given options. */
    std::vector<std::string> argv(const std::vector<std::string> & opts = {}) const {
        std::vector<std::string> args = {"vcfdist", this->query_fn, this->truth_fn, this->ref_fn};
        args.insert(args.end(), opts.begin(), opts.end());
        return args;
    }

    /** @brief Returns the path of a named entry within the temporary directory. */
    std::string path(const std::string & name) const { return this->dir.path(name); }

    /** @brief Returns the absolute path of the reference FASTA. */
    std::string ref() const { return this->ref_fn; }

private:
    TempDir dir;          ///< Owns the written VCFs and anything a case writes
    std::string query_fn; ///< Path of the written query VCF
    std::string truth_fn; ///< Path of the written truth VCF
    std::string ref_fn;   ///< Path of the checked-in reference FASTA
};

/**
 * @class InTempCwd
 * @brief Switches the working directory for the object's lifetime, creating it if needed.
 *
 * The relative-prefix cases make parse_args() create a directory beneath the current one, so they
 * have to run somewhere disposable rather than in the build tree.
 */
class InTempCwd {
public:
    /** @brief Creates the directory if absent and makes it the working directory. */
    explicit InTempCwd(const std::string & dir) : saved(std::filesystem::current_path()) {
        std::filesystem::create_directories(dir);
        std::filesystem::current_path(dir);
    }

    /** @brief Restores the previous working directory. */
    ~InTempCwd() { std::filesystem::current_path(this->saved); }

    InTempCwd(const InTempCwd &) = delete;
    InTempCwd & operator=(const InTempCwd &) = delete;

private:
    std::filesystem::path saved; ///< Working directory captured at construction
};

/**
 * @brief Calls Globals::parse_args over a writable copy of the given arguments.
 * @param[in] args Full argument vector, including argv[0]
 * @note parse_args() takes char**, so each argument is copied into its own mutable buffer. It also
 *       fopen()s the reference FASTA without ever closing it, so the descriptor is closed here;
 *       one leak per case would otherwise accumulate across the whole group.
 */
void parse(const std::vector<std::string> & args) {
    std::vector< std::vector<char> > bufs;
    for (const std::string & arg : args)
        bufs.push_back(std::vector<char>(arg.c_str(), arg.c_str() + arg.size() + 1));
    std::vector<char *> argv;
    for (std::vector<char> & buf : bufs) argv.push_back(buf.data());

    g.parse_args(static_cast<int>(argv.size()), argv.data());

    // every path that reaches the end of parse_args has opened the reference, since a failed
    // fopen() errors out instead of returning
    fclose(g.ref_fasta_fp);
    g.ref_fasta_fp = NULL;
}

/**
 * @brief Redirects stdout onto stderr and then parses, for use inside a death-test statement.
 * @param[in] args Full argument vector, including argv[0]
 * @throws ERROR if stdout cannot be redirected
 * @note A death-test matcher only sees the child's stderr, but print_usage(), print_version() and
 *       print_citation() write to stdout; without this the exit code would be all an exiting case
 *       could assert on, and every member of the argc < 4 group would look alike.
 */
void parse_showing_stdout(const std::vector<std::string> & args) {
    fflush(stdout);
    if (dup2(STDERR_FILENO, STDOUT_FILENO) < 0) ERROR("Failed to redirect stdout");
    parse(args);
}

/**
 * @brief Parses the given arguments with stderr redirected to a file, and returns what it holds.
 * @param[in] args Full argument vector, including argv[0]
 * @param[in] fn Path stderr is redirected to for the duration of the parse
 * @return Everything WARN() and INFO() wrote during the parse
 */
std::string parse_capturing_stderr(const std::vector<std::string> & args, const std::string & fn) {
    {
        StderrToFile redirect(fn);
        parse(args);
    }
    return read_text(fn);
}

/**
 * @brief Parses the given arguments with stdout captured, and returns what was captured.
 * @param[in] args Full argument vector, including argv[0]
 * @return Everything the print_* helpers wrote during the parse
 */
std::string parse_capturing_stdout(const std::vector<std::string> & args) {
    testing::internal::CaptureStdout();
    parse(args);
    return testing::internal::GetCapturedStdout();
}

/* parse_args: argc < 4 short-circuit *************************************************************/

TEST(ParseArgs, Argc1Usage) {
    GlobalsGuard guard;

    // with no arguments at all the usage message is unconditional, and the exit is a success
    EXPECT_EXIT(parse_showing_stdout({"vcfdist"}), testing::ExitedWithCode(0), "Usage: vcfdist");
}

TEST(ParseArgs, HelpShort) {
    GlobalsGuard guard;
    EXPECT_EXIT(parse_showing_stdout({"vcfdist", "-h"}), testing::ExitedWithCode(0),
            "Usage: vcfdist");
}

TEST(ParseArgs, HelpLong) {
    GlobalsGuard guard;
    EXPECT_EXIT(parse_showing_stdout({"vcfdist", "--help"}), testing::ExitedWithCode(0),
            "Usage: vcfdist");
}

TEST(ParseArgs, VersionShort) {
    GlobalsGuard guard;

    // with fewer than four arguments '-v' means version; in the main loop it means verbosity
    EXPECT_EXIT(parse_showing_stdout({"vcfdist", "-v"}), testing::ExitedWithCode(0),
            "vcfdist v" + Globals::VERSION);
}

TEST(ParseArgs, VersionLong) {
    GlobalsGuard guard;
    EXPECT_EXIT(parse_showing_stdout({"vcfdist", "--version"}), testing::ExitedWithCode(0),
            "vcfdist v" + Globals::VERSION);
}

TEST(ParseArgs, Citation) {
    GlobalsGuard guard;

    // the citation prints alone: usage only follows when help was also requested or argc is 1
    EXPECT_EXIT(parse_showing_stdout({"vcfdist", "-ci"}), testing::ExitedWithCode(0),
            "MLA Format:");
}

TEST(ParseArgs, UnknownShortArgc) {
    GlobalsGuard guard;

    // an unrecognized argument here warns and falls back to usage, but still exits successfully;
    // the same argument past the mandatory three exits 1 instead
    EXPECT_EXIT(parse_showing_stdout({"vcfdist", "--bogus"}), testing::ExitedWithCode(0),
            "Invalid usage.");
}

/* parse_args: mandatory positional arguments *****************************************************/

TEST(ParseArgs, AllMandatoryOk) {
    GlobalsGuard guard;
    ArgsFixture f;

    parse(f.argv());

    EXPECT_EQ(f.path("query.vcf"), g.query_vcf_fn);
    EXPECT_EQ(f.path("truth.vcf"), g.truth_vcf_fn);
    EXPECT_EQ(f.ref(), g.ref_fasta_fn);

    // the whole command line is recorded verbatim, argv[0] included, for the output files
    EXPECT_EQ("vcfdist " + g.query_vcf_fn + " " + g.truth_vcf_fn + " " + g.ref_fasta_fn, g.cmd);
}

TEST(ParseArgs, OptionalBeforeMandatoryWarns) {
    GlobalsGuard guard;
    ArgsFixture f;

    // the warning and usage message do not stop the parse: the flag is still taken as the query
    // filename, which then fails to open
    EXPECT_EXIT(parse_showing_stdout({"vcfdist", "-n", f.path("truth.vcf"), f.ref()}),
            testing::ExitedWithCode(1),
            "Optional arguments should be provided AFTER mandatory arguments");
}

TEST(ParseArgs, QueryOpenFailErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse({"vcfdist", f.path("absent.vcf"), f.path("truth.vcf"), f.ref()}),
            testing::ExitedWithCode(1), "Failed to open query VCF file");
}

TEST(ParseArgs, TruthOpenFailErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse({"vcfdist", f.path("query.vcf"), f.path("absent.vcf"), f.ref()}),
            testing::ExitedWithCode(1), "Failed to open truth VCF file");
}

TEST(ParseArgs, RefOpenFailErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse({"vcfdist", f.path("query.vcf"), f.path("truth.vcf"),
            f.path("absent.fasta")}), testing::ExitedWithCode(1),
            "Failed to open reference FASTA file");
}

/* parse_args: verbosity pre-pass *****************************************************************/

TEST(ParseArgs, VerbosityValid) {
    GlobalsGuard guard;
    ArgsFixture f;

    for (int v = 0; v <= 2; v++) {
        // verbosity 1 and 2 echo the command to stderr, which the redirect swallows
        StderrToFile redirect(f.path("log.txt"));
        parse(f.argv({"-v", std::to_string(v)}));
        EXPECT_EQ(v, g.verbosity);
    }
}

TEST(ParseArgs, VerbosityMissingNotReached) {
    GlobalsGuard guard;
    ArgsFixture f;

    // DOCUMENTS CURRENT BEHAVIOR, does not enforce it: the pre-pass loop stops at argc-2, so a
    // trailing '-v' is never examined and the ERROR at globals.cpp:92 cannot fire. The main loop
    // then consumes '-v' plus a value that is not there, and the parse succeeds with verbosity
    // unchanged rather than reporting the missing value.
    parse(f.argv({"-v"}));

    EXPECT_EQ(0, g.verbosity);
}

TEST(ParseArgs, VerbosityNonNumericErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-v", "high"})), testing::ExitedWithCode(1),
            "Invalid printing verbosity provided");
}

TEST(ParseArgs, VerbosityOobLowErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-v", "-1"})), testing::ExitedWithCode(1), "not a valid option");
}

TEST(ParseArgs, VerbosityOobHighErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    // the valid set is 0, 1 and 2, so 3 is one past the top
    EXPECT_EXIT(parse(f.argv({"-v", "3"})), testing::ExitedWithCode(1), "not a valid option");
}

/* parse_args: -b/--bed ***************************************************************************/

TEST(ParseArgs, BedOk) {
    GlobalsGuard guard;
    ArgsFixture f;

    parse(f.argv({"-b", data_path("tiny.bed")}));

    EXPECT_EQ(data_path("tiny.bed"), g.bed_fn);
    EXPECT_TRUE(g.bed_exists);
    ASSERT_EQ(size_t(1), g.bed.contigs.size());
    EXPECT_EQ("chr1", g.bed.contigs[0]);
    EXPECT_EQ(6, g.bed.size);
}

TEST(ParseArgs, BedMissingErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-b"})), testing::ExitedWithCode(1),
            "Option '-b' used without providing BED filename");
}

TEST(ParseArgs, BedBadFileErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    // the bedData constructor errors out on an unopenable file, so the surrounding catch block
    // never runs and the message is the constructor's
    EXPECT_EXIT(parse(f.argv({"-b", f.path("absent.bed")})), testing::ExitedWithCode(1),
            "Failed to open BED file");
}

TEST(ParseArgs, BedMalformedLineErrors) {
    GlobalsGuard guard;
    ArgsFixture f;
    const std::string bed_fn = f.path("bad.bed");
    std::ofstream(bed_fn) << "chr1\tstart\tstop\n";

    // a non-numeric coordinate is reported by the constructor too, naming the offending field and
    // line rather than the filename, so the catch block is now unreachable for every BED input
    EXPECT_EXIT(parse(f.argv({"-b", bed_fn})), testing::ExitedWithCode(1),
            "Invalid coordinate 'start' on line 1 of BED file");
}

/* parse_args: -p/--prefix ************************************************************************/

TEST(ParseArgs, PrefixRelativeDotSlash) {
    GlobalsGuard guard;
    ArgsFixture f;
    InTempCwd cwd(f.path("run"));

    // a bare relative prefix gains a "./" so that later concatenation cannot produce a path
    // starting with a flag-like character
    parse(f.argv({"-p", "out/pre"}));

    EXPECT_EQ("./out/pre", g.out_prefix);
    EXPECT_TRUE(std::filesystem::is_directory(f.path("run/out")));
}

TEST(ParseArgs, PrefixAbsoluteKept) {
    GlobalsGuard guard;
    ArgsFixture f;

    parse(f.argv({"-p", f.path("abs/pre")}));

    EXPECT_EQ(f.path("abs/pre"), g.out_prefix);
    EXPECT_TRUE(std::filesystem::is_directory(f.path("abs")));
}

TEST(ParseArgs, PrefixDotSlashKept) {
    GlobalsGuard guard;
    ArgsFixture f;
    InTempCwd cwd(f.path("run"));

    parse(f.argv({"-p", "./rel/pre"}));

    EXPECT_EQ("./rel/pre", g.out_prefix);
    EXPECT_TRUE(std::filesystem::is_directory(f.path("run/rel")));
}

TEST(ParseArgs, PrefixDotDotKept) {
    GlobalsGuard guard;
    ArgsFixture f;
    InTempCwd cwd(f.path("run/inner"));

    parse(f.argv({"-p", "../up/pre"}));

    EXPECT_EQ("../up/pre", g.out_prefix);
    EXPECT_TRUE(std::filesystem::is_directory(f.path("run/up")));
}

TEST(ParseArgs, PrefixMissingErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-p"})), testing::ExitedWithCode(1),
            "Option '-p' used without providing prefix for storing results");
}

/* parse_args: -f/--filter ************************************************************************/

TEST(ParseArgs, FilterSingle) {
    GlobalsGuard guard;
    ArgsFixture f;
    g.filters.clear();
    g.filter_ids.clear();

    parse(f.argv({"-f", "PASS"}));

    EXPECT_EQ(std::vector<std::string>({"PASS"}), g.filters);

    // one placeholder id per filter, resolved against the VCF header later
    EXPECT_EQ(std::vector<int>({-1}), g.filter_ids);
}

TEST(ParseArgs, FilterCommaList) {
    GlobalsGuard guard;
    ArgsFixture f;
    g.filters.clear();
    g.filter_ids.clear();

    parse(f.argv({"-f", "PASS,LowQual,dip"}));

    EXPECT_EQ(std::vector<std::string>({"PASS", "LowQual", "dip"}), g.filters);
    EXPECT_EQ(size_t(3), g.filter_ids.size());
}

TEST(ParseArgs, FilterTrailingComma) {
    GlobalsGuard guard;
    ArgsFixture f;
    g.filters.clear();
    g.filter_ids.clear();

    parse(f.argv({"-f", "PASS,"}));

    // the empty field after the trailing comma names no filter, so it is dropped
    EXPECT_EQ(std::vector<std::string>({"PASS"}), g.filters);

    // the id vector must stay in lockstep with the name vector
    EXPECT_EQ(std::vector<int>({-1}), g.filter_ids);
}

TEST(ParseArgs, FilterLeadingComma) {
    GlobalsGuard guard;
    ArgsFixture f;
    g.filters.clear();
    g.filter_ids.clear();

    parse(f.argv({"-f", ",PASS"}));

    EXPECT_EQ(std::vector<std::string>({"PASS"}), g.filters);
    EXPECT_EQ(std::vector<int>({-1}), g.filter_ids);
}

TEST(ParseArgs, FilterInteriorEmptyField) {
    GlobalsGuard guard;
    ArgsFixture f;
    g.filters.clear();
    g.filter_ids.clear();

    parse(f.argv({"-f", "PASS,,LowQual"}));

    EXPECT_EQ(std::vector<std::string>({"PASS", "LowQual"}), g.filters);
    EXPECT_EQ(std::vector<int>({-1, -1}), g.filter_ids);
}

TEST(ParseArgs, FilterAccumulatesAcrossFlags) {
    GlobalsGuard guard;
    ArgsFixture f;
    g.filters.clear();
    g.filter_ids.clear();

    parse(f.argv({"-f", "PASS,", "-f", "LowQual"}));

    EXPECT_EQ(std::vector<std::string>({"PASS", "LowQual"}), g.filters);
    EXPECT_EQ(std::vector<int>({-1, -1}), g.filter_ids);
}

TEST(ParseArgs, FilterMissingErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-f"})), testing::ExitedWithCode(1),
            "Option '--filter' used without providing filters");
}

TEST(ParseArgs, FilterAllEmptyFieldsErrors) {
    GlobalsGuard guard;
    ArgsFixture f;
    g.filters.clear();
    g.filter_ids.clear();

    // an empty g.filters keeps every variant, so erroring beats honoring the request backwards
    EXPECT_EXIT(parse(f.argv({"-f", ""})), testing::ExitedWithCode(1),
            "Option '--filter' provided no filter names");
}

TEST(ParseArgs, FilterOnlyCommasErrors) {
    GlobalsGuard guard;
    ArgsFixture f;
    g.filters.clear();
    g.filter_ids.clear();

    EXPECT_EXIT(parse(f.argv({"-f", ","})), testing::ExitedWithCode(1),
            "Option '--filter' provided no filter names");
}

TEST(ParseArgs, FilterEmptyAfterEarlierFlagErrors) {
    GlobalsGuard guard;
    ArgsFixture f;
    g.filters.clear();
    g.filter_ids.clear();

    // the guard is per-flag, so an earlier flag's names must not excuse a later empty one
    EXPECT_EXIT(parse(f.argv({"-f", "PASS", "-f", ""})), testing::ExitedWithCode(1),
            "Option '--filter' provided no filter names");
}

/* parse_args: -l/--largest-variant ***************************************************************/

TEST(ParseArgs, LargestOk) {
    GlobalsGuard guard;
    ArgsFixture f;

    parse(f.argv({"-l", "2000"}));

    EXPECT_EQ(2000, g.max_size);
}

TEST(ParseArgs, LargestMissingErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-l"})), testing::ExitedWithCode(1),
            "Option '--largest-variant' used without providing maximum variant size");
}

TEST(ParseArgs, LargestNonNumericErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-l", "big"})), testing::ExitedWithCode(1),
            "Invalid maximum variant size provided");
}

/* parse_args: -sv/--sv-threshold *****************************************************************/

TEST(ParseArgs, SvOk) {
    GlobalsGuard guard;
    ArgsFixture f;

    parse(f.argv({"-sv", "100"}));

    EXPECT_EQ(100, g.sv_threshold);
}

TEST(ParseArgs, SvMissingErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-sv"})), testing::ExitedWithCode(1),
            "Option '--sv-threshold' used without providing an SV threshold size");
}

TEST(ParseArgs, SvNonNumericErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-sv", "large"})), testing::ExitedWithCode(1),
            "Invalid SV threshold size provided");
}

TEST(ParseArgs, SvTooSmallErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    // a threshold of 1 would classify every SNP as an SV, so 2 is the smallest accepted
    EXPECT_EXIT(parse(f.argv({"-sv", "1"})), testing::ExitedWithCode(1),
            "Must provide larger SV threshold size");
}

/* parse_args: -q/--min-qual **********************************************************************/

TEST(ParseArgs, MinQualOk) {
    GlobalsGuard guard;
    ArgsFixture f;

    parse(f.argv({"-q", "10"}));

    EXPECT_EQ(10, g.min_qual);
}

TEST(ParseArgs, MinQualMissingErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-q"})), testing::ExitedWithCode(1),
            "Option '--min-qual' used without providing minimum variant quality");
}

TEST(ParseArgs, MinQualNonNumericErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-q", "high"})), testing::ExitedWithCode(1),
            "Invalid minimum variant quality provided");
}

TEST(ParseArgs, MinQualNegativeErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-q", "-1"})), testing::ExitedWithCode(1),
            "Must provide non-negative minimum variant quality");
}

/* parse_args: -mq/--max-qual *********************************************************************/

TEST(ParseArgs, MaxQualOk) {
    GlobalsGuard guard;
    ArgsFixture f;

    parse(f.argv({"-mq", "99"}));

    EXPECT_EQ(99, g.max_qual);
}

TEST(ParseArgs, MaxQualMissingErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-mq"})), testing::ExitedWithCode(1),
            "Option '--max-qual' used without providing maximum variant quality");
}

TEST(ParseArgs, MaxQualNonNumericErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-mq", "high"})), testing::ExitedWithCode(1),
            "Invalid maximum variant quality provided");
}

/* parse_args: -n/--no-output-files ***************************************************************/

TEST(ParseArgs, NoOutput) {
    GlobalsGuard guard;
    ArgsFixture f;
    ASSERT_TRUE(g.write);

    // the flag takes no value, so the loop advances by one rather than two
    parse(f.argv({"-n"}));

    EXPECT_FALSE(g.write);
}

/* parse_args: -x/--mismatch-penalty **************************************************************/

TEST(ParseArgs, MismatchOk) {
    GlobalsGuard guard;
    ArgsFixture f;

    parse(f.argv({"-x", "7"}));

    EXPECT_EQ(7, g.sub);
}

TEST(ParseArgs, MismatchMissingErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-x"})), testing::ExitedWithCode(1),
            "Option '-x' used without providing mismatch penalty");
}

TEST(ParseArgs, MismatchNegativeErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-x", "-1"})), testing::ExitedWithCode(1),
            "Must provide non-negative mismatch penalty");
}

/* parse_args: -o/--gap-open-penalty **************************************************************/

TEST(ParseArgs, GapOpenOk) {
    GlobalsGuard guard;
    ArgsFixture f;

    parse(f.argv({"-o", "9"}));

    EXPECT_EQ(9, g.open);
}

TEST(ParseArgs, GapOpenMissingErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-o"})), testing::ExitedWithCode(1),
            "Option '-o' used without providing gap-opening penalty");
}

TEST(ParseArgs, GapOpenNegativeErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-o", "-1"})), testing::ExitedWithCode(1),
            "Must provide non-negative gap-opening penalty");
}

/* parse_args: -e/--gap-extend-penalty ************************************************************/

TEST(ParseArgs, GapExtendOk) {
    GlobalsGuard guard;
    ArgsFixture f;

    parse(f.argv({"-e", "3"}));

    EXPECT_EQ(3, g.extend);
}

TEST(ParseArgs, GapExtendMissingErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-e"})), testing::ExitedWithCode(1),
            "Option '-e' used without providing gap-extension penalty");
}

TEST(ParseArgs, GapExtendNegativeErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-e", "-1"})), testing::ExitedWithCode(1),
            "Must provide non-negative gap-extension penalty");
}

/* parse_args: -i/--max-iterations ****************************************************************/

TEST(ParseArgs, IterationsOk) {
    GlobalsGuard guard;
    ArgsFixture f;

    // issue #46 removes this flag; delete this case along with it
    parse(f.argv({"-i", "3"}));

    EXPECT_EQ(3, g.max_cluster_itrs);
}

TEST(ParseArgs, IterationsTooSmallErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    // issue #46 removes this flag; delete this case along with it
    EXPECT_EXIT(parse(f.argv({"-i", "0"})), testing::ExitedWithCode(1),
            "Max cluster iterations must be positive");
}

/* parse_args: -s/--max-supercluster-size *********************************************************/

TEST(ParseArgs, SuperclusterOk) {
    GlobalsGuard guard;
    ArgsFixture f;

    // issue #47 renames this flag to -sc; retarget this case along with it
    parse(f.argv({"-s", "20000"}));

    EXPECT_EQ(20000, g.max_supercluster_size);
}

TEST(ParseArgs, SuperclusterTooSmallErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    // issue #47 renames this flag to -sc; retarget this case along with it
    EXPECT_EXIT(parse(f.argv({"-s", "0"})), testing::ExitedWithCode(1),
            "Max supercluster size must be positive");
}

/* parse_args: -t/--max-threads *******************************************************************/

TEST(ParseArgs, ThreadsOk) {
    GlobalsGuard guard;
    ArgsFixture f;

    parse(f.argv({"-t", "4"}));

    EXPECT_EQ(4, g.max_threads);

    // the scheduler is recomputed at the end of the parse, so the steps follow the new count
    EXPECT_EQ(std::vector<int>({4, 2, 1}), g.thread_steps);
    EXPECT_EQ(3, g.thread_nsteps);
}

TEST(ParseArgs, ThreadsMissingErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-t"})), testing::ExitedWithCode(1),
            "Option '--max-threads' used without providing max threads");
}

TEST(ParseArgs, ThreadsNonNumericErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-t", "many"})), testing::ExitedWithCode(1),
            "Invalid max threads provided");
}

TEST(ParseArgs, ThreadsTooSmallErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    // zero threads would leave the scheduler with no steps at all
    EXPECT_EXIT(parse(f.argv({"-t", "0"})), testing::ExitedWithCode(1),
            "Max threads must be positive");
}

/* parse_args: -ct/--credit-threshold *************************************************************/

TEST(ParseArgs, CreditOk) {
    GlobalsGuard guard;
    ArgsFixture f;

    parse(f.argv({"-ct", "0.5"}));

    EXPECT_DOUBLE_EQ(0.5, g.credit_threshold);
}

TEST(ParseArgs, CreditMissingErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-ct"})), testing::ExitedWithCode(1),
            "Option '--credit-threshold' used without providing value");
}

TEST(ParseArgs, CreditNonNumericErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-ct", "most"})), testing::ExitedWithCode(1),
            "Invalid credit threshold provided");
}

TEST(ParseArgs, CreditZeroErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    // the interval is open at the bottom, so zero is rejected rather than meaning "no threshold"
    EXPECT_EXIT(parse(f.argv({"-ct", "0"})), testing::ExitedWithCode(1), "must be on the interval");
}

TEST(ParseArgs, CreditAboveOneErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-ct", "1.5"})), testing::ExitedWithCode(1),
            "must be on the interval");
}

TEST(ParseArgs, CreditUpperInclusive) {
    GlobalsGuard guard;
    ArgsFixture f;

    // the interval is closed at the top, so exact credit is a valid requirement
    parse(f.argv({"-ct", "1.0"}));

    EXPECT_DOUBLE_EQ(1.0, g.credit_threshold);
}

/* parse_args: -r/--max-ram ***********************************************************************/

TEST(ParseArgs, RamOk) {
    GlobalsGuard guard;
    ArgsFixture f;

    parse(f.argv({"-r", "16.5"}));

    EXPECT_DOUBLE_EQ(16.5, g.max_ram);
}

TEST(ParseArgs, RamMissingErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-r"})), testing::ExitedWithCode(1),
            "Option '--max-ram' used without providing max RAM");
}

TEST(ParseArgs, RamNonNumericErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-r", "lots"})), testing::ExitedWithCode(1),
            "Invalid max RAM provided");
}

TEST(ParseArgs, RamTrailingUnitsWarns) {
    GlobalsGuard guard;
    ArgsFixture f;

    // the units are implicitly GB, so a written unit is a warning rather than an error and the
    // leading number is still used
    const std::string log = parse_capturing_stderr(f.argv({"-r", "64GB"}), f.path("log.txt"));

    EXPECT_DOUBLE_EQ(64.0, g.max_ram);
    EXPECT_TRUE(logged(log, "Trailing text 'GB' ignored for --max-ram")) << log;
}

TEST(ParseArgs, RamZeroErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"-r", "0"})), testing::ExitedWithCode(1), "Max RAM must be positive");
}

TEST(ParseArgs, RamNegativeErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    // zero and negative share one guard, since neither leaves any RAM to divide among threads
    EXPECT_EXIT(parse(f.argv({"-r", "-8"})), testing::ExitedWithCode(1),
            "Max RAM must be positive");
}

/* parse_args: main-loop help, version, citation, and unknown options *****************************/

TEST(ParseArgs, HelpMainLoop) {
    GlobalsGuard guard;
    ArgsFixture f;

    // past the mandatory arguments the usage message is deferred to the end of the parse, and the
    // parse returns instead of exiting
    const std::string out = parse_capturing_stdout(f.argv({"-h"}));

    EXPECT_NE(std::string::npos, out.find("Usage: vcfdist"));
}

TEST(ParseArgs, VersionMainLoop) {
    GlobalsGuard guard;
    ArgsFixture f;

    // only the long form is a version request here; bare '-v' is verbosity
    const std::string out = parse_capturing_stdout(f.argv({"--version"}));

    EXPECT_NE(std::string::npos, out.find("vcfdist v" + Globals::VERSION));
}

TEST(ParseArgs, CitationMainLoop) {
    GlobalsGuard guard;
    ArgsFixture f;

    const std::string out = parse_capturing_stdout(f.argv({"-ci"}));

    EXPECT_NE(std::string::npos, out.find("MLA Format:"));
    EXPECT_NE(std::string::npos, out.find("BibTeX Format:"));
}

TEST(ParseArgs, UnknownOptionErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    EXPECT_EXIT(parse(f.argv({"--bogus"})), testing::ExitedWithCode(1), "Unexpected option");
}

TEST(ParseArgs, VerbositySkippedInMainLoop) {
    GlobalsGuard guard;
    ArgsFixture f;

    // the pre-pass already applied the value, so the main loop consumes both tokens without effect
    // and keeps parsing what follows
    const std::string log = parse_capturing_stderr(f.argv({"-v", "2", "-n"}), f.path("log.txt"));

    EXPECT_EQ(2, g.verbosity);
    EXPECT_FALSE(g.write);
    EXPECT_TRUE(logged(log, "Command:")) << log;
}

/* parse_args: cross-field validation *************************************************************/

TEST(ParseArgs, MaxQualLtMinQualErrors) {
    GlobalsGuard guard;
    ArgsFixture f;

    // checked after the loop, so the order the two flags appear in does not matter
    EXPECT_EXIT(parse(f.argv({"-q", "50", "-mq", "10"})), testing::ExitedWithCode(1),
            "Maximum variant quality must exceed minimum variant quality");
}

TEST(ParseArgs, MaxSizeLtSvWarns) {
    GlobalsGuard guard;
    ArgsFixture f;

    // every variant small enough to evaluate is below the SV threshold, so the SV category is empty
    const std::string log = parse_capturing_stderr(f.argv({"-l", "20"}), f.path("log.txt"));

    EXPECT_TRUE(logged(log, "No SVs will be evaluated")) << log;
}

TEST(ParseArgs, SuperclusterLtMaxSizePlus2Errors) {
    GlobalsGuard guard;
    ArgsFixture f;

    // a supercluster has to hold the largest evaluated variant plus one flanking base on each side
    EXPECT_EXIT(parse(f.argv({"-l", "200", "-s", "100"})), testing::ExitedWithCode(1),
            "Invalid option selected: --max-supercluster-size");
}

TEST(ParseArgs, MaxSize1WarnsSnpsOnly) {
    GlobalsGuard guard;
    ArgsFixture f;

    const std::string log = parse_capturing_stderr(f.argv({"-l", "1"}), f.path("log.txt"));

    // a maximum size of 1 is below every valid SV threshold, so both warnings necessarily fire
    EXPECT_TRUE(logged(log, "Only SNPs will be evaluated")) << log;
    EXPECT_TRUE(logged(log, "No SVs will be evaluated")) << log;
}

/* set_thread_ram_steps ***************************************************************************/

TEST(SetThreadRamSteps, Default64) {
    GlobalsGuard guard;
    g.max_threads = 64;

    g.set_thread_ram_steps();

    EXPECT_EQ(std::vector<int>({64, 32, 16, 8, 4, 2, 1}), g.thread_steps);
    EXPECT_EQ(7, g.thread_nsteps);
}

TEST(SetThreadRamSteps, SingleThread) {
    GlobalsGuard guard;
    g.max_threads = 1;

    g.set_thread_ram_steps();

    // halving 1 yields 0, which ends the loop, so a single-threaded run has one step
    EXPECT_EQ(std::vector<int>({1}), g.thread_steps);
    EXPECT_EQ(1, g.thread_nsteps);
}

TEST(SetThreadRamSteps, NonPowerOfTwo) {
    GlobalsGuard guard;
    g.max_threads = 10;

    g.set_thread_ram_steps();

    // the halving is integer division, so 5 becomes 2 rather than 2.5
    EXPECT_EQ(std::vector<int>({10, 5, 2, 1}), g.thread_steps);
    EXPECT_EQ(4, g.thread_nsteps);
}

TEST(SetThreadRamSteps, ZeroThreadsEmpty) {
    GlobalsGuard guard;
    g.max_threads = 0;

    // parse_args rejects a thread count below 1, but the loop guard is what makes an out-of-range
    // count yield no steps rather than dividing by zero
    g.set_thread_ram_steps();

    EXPECT_TRUE(g.thread_steps.empty());
    EXPECT_TRUE(g.ram_steps.empty());
    EXPECT_EQ(0, g.thread_nsteps);
}

TEST(SetThreadRamSteps, RamDividedPerThread) {
    GlobalsGuard guard;
    g.max_threads = 8;
    g.max_ram = 16;

    g.set_thread_ram_steps();

    // fewer threads each get a larger share of the same total, which is the point of the steps
    ASSERT_EQ(g.thread_steps.size(), g.ram_steps.size());
    for (size_t i = 0; i < g.ram_steps.size(); i++) {
        EXPECT_FLOAT_EQ(float(g.max_ram / g.thread_steps[i]), g.ram_steps[i]) << "step " << i;
    }
    EXPECT_EQ(std::vector<float>({2, 4, 8, 16}), g.ram_steps);
}

TEST(SetThreadRamSteps, RepeatedCallsReplace) {
    GlobalsGuard guard;
    g.max_threads = 8;

    // both vectors are cleared first, so the constructor's steps and the parse's do not accumulate
    g.set_thread_ram_steps();
    const std::vector<int> first = g.thread_steps;
    g.set_thread_ram_steps();

    EXPECT_EQ(first, g.thread_steps);
    EXPECT_EQ(4, g.thread_nsteps);
}

TEST(SetThreadRamSteps, NstepsMatchesVectorSize) {
    GlobalsGuard guard;

    // the count is maintained alongside the vectors, and callers index the vectors with it
    for (int threads : {1, 3, 8, 100}) {
        g.max_threads = threads;
        g.set_thread_ram_steps();
        EXPECT_EQ(size_t(g.thread_nsteps), g.thread_steps.size()) << threads << " threads";
        EXPECT_EQ(size_t(g.thread_nsteps), g.ram_steps.size()) << threads << " threads";
    }
}

/* print_version, print_usage, print_citation *****************************************************/

/** @brief Returns everything print_usage() writes to stdout. */
std::string usage_text() {
    testing::internal::CaptureStdout();
    g.print_usage();
    return testing::internal::GetCapturedStdout();
}

TEST(PrintVersion, Format) {
    testing::internal::CaptureStdout();

    g.print_version();

    EXPECT_EQ(Globals::PROGRAM + " v" + Globals::VERSION + "\n",
            testing::internal::GetCapturedStdout());
}

TEST(PrintUsage, RequiredSection) {
    GlobalsGuard guard;

    const std::string usage = usage_text();

    EXPECT_NE(std::string::npos,
            usage.find("Usage: vcfdist <query.vcf> <truth.vcf> <ref.fasta> [options]"));
    EXPECT_NE(std::string::npos, usage.find("Required:"));
    EXPECT_NE(std::string::npos, usage.find("query.vcf"));
    EXPECT_NE(std::string::npos, usage.find("truth.vcf"));
    EXPECT_NE(std::string::npos, usage.find("ref.fasta"));
}

TEST(PrintUsage, ListsDocumentedFlags) {
    GlobalsGuard guard;

    const std::string usage = usage_text();

    for (const std::string & flag : {"-b, --bed", "-v, --verbosity", "-p, --prefix",
            "-n, --no-output-files", "-f, --filter", "-l, --largest-variant",
            "-sv, --sv-threshold", "-q, --min-qual", "-mq, --max-qual",
            "-s, --max-supercluster-size", "-ct, --credit-threshold", "-t, --max-threads",
            "-r, --max-ram", "-h, --help", "-ci, --citation", "-v, --version"}) {
        EXPECT_NE(std::string::npos, usage.find(flag)) << "undocumented flag '" << flag << "'";
    }
}

TEST(PrintUsage, OmitsHiddenFlags) {
    GlobalsGuard guard;

    const std::string usage = usage_text();

    // these four are parsed but their usage lines are commented out, so they are tuning knobs
    // rather than part of the advertised interface
    for (const std::string & flag : {"--max-iterations", "--mismatch-penalty",
            "--gap-open-penalty", "--gap-extend-penalty"}) {
        EXPECT_EQ(std::string::npos, usage.find(flag)) << "advertised hidden flag '" << flag << "'";
    }
}

TEST(PrintUsage, InterpolatesDefaults) {
    GlobalsGuard guard;
    g.verbosity = 2;
    g.max_size = 123;
    g.sv_threshold = 45;
    g.min_qual = 7;
    g.max_qual = 88;
    g.max_supercluster_size = 999;
    g.credit_threshold = 0.5;
    g.max_threads = 3;
    g.max_ram = 12.5;

    // the bracketed values are the live settings, not compiled-in literals
    const std::string usage = usage_text();

    EXPECT_NE(std::string::npos, usage.find("--verbosity <INTEGER> [2]"));
    EXPECT_NE(std::string::npos, usage.find("--largest-variant <INTEGER> [123]"));
    EXPECT_NE(std::string::npos, usage.find("--sv-threshold <INTEGER> [45]"));
    EXPECT_NE(std::string::npos, usage.find("--min-qual <INTEGER> [7]"));
    EXPECT_NE(std::string::npos, usage.find("--max-qual <INTEGER> [88]"));
    EXPECT_NE(std::string::npos, usage.find("--max-supercluster-size <INTEGER> [999]"));
    EXPECT_NE(std::string::npos, usage.find("--credit-threshold <FLOAT> [0.50]"));
    EXPECT_NE(std::string::npos, usage.find("--max-threads <INTEGER> [3]"));
    EXPECT_NE(std::string::npos, usage.find("--max-ram <FLOAT> [12.50GB]"));
}

TEST(PrintCitation, BothFormats) {
    testing::internal::CaptureStdout();

    g.print_citation();
    const std::string cite = testing::internal::GetCapturedStdout();

    EXPECT_NE(std::string::npos, cite.find("MLA Format:"));
    EXPECT_NE(std::string::npos, cite.find("BibTeX Format:"));
    EXPECT_NE(std::string::npos, cite.find("@article{dunn2024vcfdist,"));
    EXPECT_NE(std::string::npos, cite.find("Genome Biology"));
    EXPECT_NE(std::string::npos, cite.find("10.1186/s13059-024-03394-5"));
}

} // namespace
