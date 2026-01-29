#include "../src/core/logging.h"
#include "../src/core/timer.h"
#include "../src/core/config.h"

#include <gtest/gtest.h>
#include <sstream>
#include <thread>

using namespace fem;

// Helper to capture std::cout during a scope
struct CoutCapture {
    std::streambuf* old = nullptr;
    std::ostringstream oss;
    CoutCapture() { old = std::cout.rdbuf(oss.rdbuf()); }
    ~CoutCapture() { std::cout.rdbuf(old); }
    std::string str() const { return oss.str(); }
};

// 1) logging indentation up/down and get_indentation()
TEST(CoreLogging, IndentationUpDown) {
    // Reset indentation
    while (!fem::logging::get_indentation().empty()) fem::logging::down();
    EXPECT_EQ(fem::logging::get_indentation(), "");
    fem::logging::up();
    EXPECT_EQ(fem::logging::get_indentation(), "  ");
    fem::logging::up();
    EXPECT_EQ(fem::logging::get_indentation(), std::string(4, ' '));
    fem::logging::down();
    EXPECT_EQ(fem::logging::get_indentation(), std::string(2, ' '));
}

// 2) logging wrapping preserves leading whitespace of segment
TEST(CoreLogging, WrappingPreservesLeadingWhitespace) {
    // Ensure deterministic width
    fem::logging::set_console_width(40);
    // Reset indentation
    while (!fem::logging::get_indentation().empty()) fem::logging::down();

    // Message with leading spaces in the segment
    const std::string payload = "    This is a long line that should wrap nicely and keep the indent.";

    CoutCapture cap;
    fem::logging::info(true, payload);
    const std::string out = cap.str();

    // Expect at least 2 physical lines
    std::vector<std::string> lines;
    std::stringstream ss(out);
    std::string line;
    while (std::getline(ss, line)) if (!line.empty()) lines.push_back(line);
    ASSERT_GE(lines.size(), 2u);

    // Base prefix is "[INFO] " with zero indentation
    const std::string prefix = "[INFO] ";
    ASSERT_TRUE(lines[0].rfind(prefix, 0) == 0);
    ASSERT_TRUE(lines[1].rfind(prefix, 0) == 0);

    // After the prefix, both lines should start with the same 4 leading spaces
    EXPECT_TRUE(lines[0].substr(prefix.size()).rfind("    ", 0) == 0);
    EXPECT_TRUE(lines[1].substr(prefix.size()).rfind("    ", 0) == 0);
}

// 3) timer positive elapsed and Timer::measure result passthrough
TEST(CoreTimer, ElapsedPositiveAndMeasure) {
    Timer t{};
    t.start();
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
    t.stop();
    EXPECT_GT(t.elapsed(), 0u);

    auto res = Timer::measure([](){
        long long s = 0; for (int i = 0; i < 1000; ++i) s += i; return s;
    }, "sum");
    EXPECT_EQ(res, 999LL * 1000 / 2);
}

// 4) config defaults and mutability
TEST(CoreConfig, DefaultsAndMutability) {
    // default
    EXPECT_EQ(global_config.max_threads, 1);
    auto old = global_config.max_threads;
    global_config.max_threads = 3;
    EXPECT_EQ(global_config.max_threads, 3);
    // restore
    global_config.max_threads = old;
}

