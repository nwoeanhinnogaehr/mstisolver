#pragma once

#include <mutex>

static std::mutex _output_mutex;
static bool _log = true;

#define LOG(...) if(_log){std::lock_guard<std::mutex> _guard(_output_mutex); cerr << __VA_ARGS__; cerr << flush;}
#define OUT(...) {std::lock_guard<std::mutex> _guard(_output_mutex); cout << __VA_ARGS__; cout << flush;}
