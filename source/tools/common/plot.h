/*
MIT License

Copyright (c) 2017  Joe Hood

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once

#include <string>
#include <string_view>
#include <vector>

class AsciiPlotter {
    private:
    static constexpr size_t          MAX_CURVES = 10;
    size_t                           margin     = 10;
    size_t                           width      = 100;
    size_t                           height     = 50;
    std::string                      title      = "";
    std::string                      xlabel     = "";
    std::string                      ylabel     = "";
    bool                             legend     = false;
    std::vector<double>              xdata;
    std::vector<std::vector<double>> ydata;
    std::vector<std::string>         labels;
    std::vector<char>                markers;

    inline static std::vector<double>              xdata_s   = {};
    inline static std::vector<std::vector<double>> ydata_s   = {};
    inline static std::vector<std::string>         labels_s  = {};
    inline static std::vector<char>                markers_s = {};

    public:
    AsciiPlotter() = default;
    AsciiPlotter(std::string_view title);
    AsciiPlotter(std::string_view title, size_t width, size_t height);
    void addPlot(const std::vector<double> &xdata, const std::vector<double> &ydata, std::string_view label = "", char marker = '.');
    void addPlot(const std::vector<double> &ydata, std::string_view label = "", char marker = '.');
    void addStaticPlot(const std::vector<double> &xdata, const std::vector<double> &ydata, std::string_view label = "", char marker = '.');
    void addStaticPlot(const std::vector<double> &ydata, std::string_view label = "", char marker = '.');

    void show();
    void set_margin(size_t margin);
    void set_xlabel(std::string_view label);
    void set_ylabel(std::string_view label);
    void enable_legend();
};
