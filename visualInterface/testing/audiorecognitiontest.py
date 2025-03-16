import sounddevice as sd
import numpy as np
import json
from vosk import Model, KaldiRecognizer

SAMPLE_RATE = 16000
DURATION = 5  # seconds

# Initialize model (adjust path as needed)
model = Model("/Users/hililbby/Library/Mobile Documents/com~apple~CloudDocs/UT Austin/JM rotation/distractor_NFT/visualInterface/vosk/vosk-model-en-us-0.22")
recognizer = KaldiRecognizer(model, SAMPLE_RATE)

print("Recording for", DURATION, "seconds...")

# Record audio data
audio_data = sd.rec(int(DURATION * SAMPLE_RATE), samplerate=SAMPLE_RATE, channels=1, dtype='int16')
sd.wait()  # Wait until recording is finished

# Convert the NumPy array to bytes
audio_bytes = audio_data.tobytes()

if recognizer.AcceptWaveform(audio_bytes):
    result = recognizer.Result()
    text = json.loads(result)["text"]
    print("Final recognized text:", text)
else:
    print("No complete result, try again.")


