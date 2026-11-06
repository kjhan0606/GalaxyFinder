#!/usr/bin/env python3
"""Use the same host-unit summaries for the declared stronger-noise extension."""
import measure_image_noise_boundary as extension
extension.configure()
import summarize_image_sensitivity as summary

if __name__=="__main__":summary.main()
