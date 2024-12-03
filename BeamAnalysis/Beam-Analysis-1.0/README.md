# Beam Analysis 

This code was written by several wcte shifters over the last few weeks taking 8 hour shifts writing in the same files.
At the very least, I know [@acraplet](https://github.com/Acraplet) made the first working version using the tuples made by [@arturof's](https://github.com/arturof) code, and [@bferrazzi](https://github.com/bferrazzi) later added the python viewer.
I wanted to keep this version controlled so we can track changes, make edits without worrying about breaking things, and don't have to risk losing any of the work people have put in. 

# Usage

Run 
```
. analyse_tuple.sh $RUN_NO
```
whle replacing `$RUN_NO` with the run number. No leading zero is needed.

You will be prompted for a password; this is the normal SK password. 

# Buffering 

The buffering aspect of this is still a work in progress!

The data are first buffered by `analyse_tuple`'s main function. This sets them aside but does not add them to the historam.

Eventually, this will be setup so that data are being buffered. Then, only on proton events and where we can know with certainty that there is a match of QDC and TDC counts, the buffered events would be "flushed" into a histogram. 
When a mismatch is found, the buffer would instead then be cleared.

# Modifying 

## Compiling

You should only have to run `make` in the `beam_analysis` folder to compile this. 

## Histograming 

Filling of histograms is done in the `buffer.cpp` file by the `EventBuffer::Flush` method. 
To make a new histogram, define the new `THxD` in the `buffer.h` header, initialize it in the `EventBuffer::EventBuffer` constructor, and fill it in the `Flush` method. 
You should then draw it in the `EventBuffer::MakePlots` method. 

## Collecting more information

If you add in a new calculation to `analyse_tuple.cpp` and want to incorporate it here.
- Prepare a new `THxD` as above.
- Set up a new buffer type in the `BufferT` enum in `buffer.h`
- Add a new buffer in the `EventBuffer` class private attributes. 
- Add a new entry to the switch in the `EventBuffer::Buffer` method. If this is a per-PMT calculation, add it to the method taking a `pmt_id`. 
- Update the `Flush` method to flush your buffer into a histogram.
- Update the `Clear` method to properly clear out your buffer. 
- Update the `MakePlots` method with a canvas for your new histogram.

## Previous Versions

The previous version of this code is backed up in both `~/home/backup_beam_analysis` and also under the first commit of this repo. 