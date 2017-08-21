#example of how to run plot_beams.py
python /data/code/git/mwa_beam_measurement/plot_beams.py --generate_pb_map --plot_pb_map --obs_start_time="1501833532" \
   --obs_end_time="1501833542" --raw_data_dir="/data/code/git/mwa_beam_measurement/test_data" --ref_tile_list="rf1YY" --AUT_tile_list="052YY" \
      --ref_signal_threshold="-150" --AUT_signal_threshold="-150" --alt_threshold="10" --plot_1D
