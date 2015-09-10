extern "C" {
#include "nsgt.h"
#include "examples/example_nsgt_error.h"
#include "examples/example_nsgt_stream_error.h"
#include "examples/example_spectrogram.h"
}

#include <sndfile.h>
#include <complex>
#include <iostream>
#include <vector>

bool read_mono_wav(std::string filename, SF_INFO *sfinfo, std::vector<double>& input) {
    SNDFILE *f = sf_open(filename.c_str(), SFM_READ, sfinfo);
    if (!f || sfinfo->channels <= 0 || sfinfo->frames <= 0 || sfinfo->samplerate <= 0) {
        std::cout << "Unable to open file 'in.wav'" << std::endl;
        std::cout << "Please make sure you start the application with the 'examples' working directory" << std::endl;
		return false;
    }
    if (sfinfo->channels != 1) {
        std::cout << "Only the first channel of wav file will be used in the transform" << std::endl;
    }
    std::size_t len = (size_t) sfinfo->frames;
    input.resize(len);
    int err = 0;
    std::vector<double> buf;
    buf.resize((std::size_t) sfinfo->channels);
    for (size_t i = 0; i < len; i++) {
        int read = (int) sf_read_double(f, buf.data(), sfinfo->channels);
        if (read != sfinfo->channels) {
            err = 1;
            break;
        }
        input[i] = buf[0];
    }
    sf_close(f);
    if (err) {
        printf("Error while reading the file.\r\n");
		return false;
	}
    return true;
}

int main() {
    SF_INFO info = {};
	std::vector<double> input;

	if (read_mono_wav("in.wav", &info, input)) {
		example_nsgt_error(input.data(), input.size(), info.samplerate);
		example_stream_error(input.data(), input.size(), info.samplerate);
		example_spectrogram(input.data(), input.size(), info.samplerate);
	}

#ifdef _MSC_VER
	std::cout << "done." << std::endl;
    std::cin.ignore(); //windows closes the terminal window after the application exits
#endif

    return 0;
}
