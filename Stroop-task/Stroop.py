#!/usr/bin/env python3
"""
Stroop Experiment re-implemented using Pygame

Requirements:
    - pygame
    - pandas (if you wish to load conditions from an Excel file)
    - openpyxl (if you use an Excel file for conditions)

Instructions:
    - The experiment displays instructions on screen.
    - Press any key to start.
    - On each trial, after 0.5 sec a word appears in a given color.
    - Respond by pressing:
          left  -> red
          right -> blue
          down  -> green
    - Your reaction time and whether your response was correct are recorded.
    - Data are saved in a CSV file under the “data” folder.
"""

import pygame
import sys
import random
import time
import datetime
import os
import csv

# Try to import pandas to read conditions from an Excel file.
# If not available, the script will fall back to a default condition list.
try:
    import pandas as pd
    HAVE_PANDAS = True
except ImportError:
    HAVE_PANDAS = False

def get_color(color_str):
    """Map a color name to an RGB tuple."""
    if isinstance(color_str, str):
        lower = color_str.lower()
        if lower == "red":
            return (255, 0, 0)
        elif lower == "blue":
            return (0, 0, 255)
        elif lower == "green":
            return (0, 255, 0)
    return (255, 255, 255)  # default white

def load_conditions(filename="conditions.xlsx"):
    """
    Load conditions from an Excel file.
    The file should contain columns named (at least): word, colour, corrAns
    """
    if HAVE_PANDAS:
        try:
            df = pd.read_excel(filename)
            conditions = df.to_dict(orient="records")
            return conditions
        except Exception as e:
            print("Error reading conditions file:", e)
    # Fall back to a default list if reading fails.
    print("Using default conditions.")
    return [
        {"word": "RED",   "colour": "blue",  "corrAns": "right"},
        {"word": "BLUE",  "colour": "red",   "corrAns": "left"},
        {"word": "GREEN", "colour": "red",   "corrAns": "left"},
        {"word": "RED",   "colour": "green", "corrAns": "down"},
    ]

def main():
    # ================================
    # Get participant and session info
    # ================================
    participant = input("Enter participant ID: ")
    session = input("Enter session number: ")
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    expName = "stroopex"
    
    # Create data directory if needed and a filename for the CSV output.
    data_dir = "data"
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    data_filename = os.path.join(data_dir, f"{participant}_{expName}_{timestamp}.csv")
    
    # =====================================
    # Load conditions (and replicate 5 times)
    # =====================================
    conditions = load_conditions("conditions.xlsx")
    # replicate the conditions 5 times and randomize order
    trials = conditions * 5
    random.shuffle(trials)
    
    # ================================
    # Initialize Pygame and create window
    # ================================
    pygame.init()
    screen_width, screen_height = 1366, 768
    # Use FULLSCREEN mode; if you prefer windowed mode, remove pygame.FULLSCREEN.
    screen = pygame.display.set_mode((screen_width, screen_height), pygame.FULLSCREEN)
    pygame.display.set_caption("Stroop Experiment")
    clock = pygame.time.Clock()
    
    # Set up fonts (you can adjust the size as needed)
    font = pygame.font.SysFont("Arial", 50)
    instruction_font = pygame.font.SysFont("Arial", 40)
    
    # ======================
    # Display Instructions
    # ======================
    instructions_text = (
        "Remember you choose the color of the letters, ignoring the word:\n\n"
        "left -> red\n"
        "right -> blue\n"
        "down -> green\n\n"
        "Press any key to start."
    )
    # Split into lines for multi-line display.
    instructions_lines = instructions_text.split("\n")
    
    def draw_instructions():
        screen.fill((0, 0, 0))
        # Start drawing a little above the vertical center.
        start_y = screen_height / 2 - (len(instructions_lines) * 30)
        for i, line in enumerate(instructions_lines):
            text_surface = instruction_font.render(line, True, (255, 255, 255))
            text_rect = text_surface.get_rect(center=(screen_width / 2, start_y + i * 50))
            screen.blit(text_surface, text_rect)
        pygame.display.flip()
    
    draw_instructions()
    waiting = True
    while waiting:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_ESCAPE:
                    pygame.quit()
                    sys.exit()
                else:
                    waiting = False
        clock.tick(60)
    
    # ==========================================
    # Main Trial Loop
    # ==========================================
    # List to collect trial data.
    data_list = []
    trial_number = 0
    
    for trial in trials:
        trial_number += 1
        # Clear event queue before starting a trial.
        pygame.event.clear()
        
        # Timing: we want a 0.5 sec delay before showing the target.
        trial_start = pygame.time.get_ticks()  # in milliseconds
        target_displayed = False  # flag to know when we first show the target
        response = None
        rt = None  # reaction time
        
        # Pre-render variables (will be set when target is drawn)
        target_surface = None
        target_rect = None
        
        # Trial loop: run until a response is recorded or the trial times out.
        trial_done = False
        while not trial_done:
            current_time = pygame.time.get_ticks()
            elapsed = (current_time - trial_start) / 1000.0  # elapsed time in seconds
            
            # Process events.
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    pygame.quit()
                    sys.exit()
                elif event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_ESCAPE:
                        pygame.quit()
                        sys.exit()
                    # Only record response if the target is on screen.
                    if target_displayed and response is None:
                        if event.key in (pygame.K_LEFT, pygame.K_RIGHT, pygame.K_DOWN):
                            response = event.key
                            # Reaction time measured from when target appeared (after 0.5 sec).
                            rt = elapsed - 0.5
                            trial_done = True
            
            # Draw the trial screen.
            screen.fill((0, 0, 0))
            if elapsed >= 0.5:
                if not target_displayed:
                    # Once 0.5 sec have passed, render the target stimulus.
                    # The trial dictionary should contain "word" and "colour".
                    color = get_color(trial["colour"])
                    target_surface = font.render(trial["word"], True, color)
                    target_rect = target_surface.get_rect(center=(screen_width / 2, screen_height / 2))
                    target_displayed = True
                # Display the target for up to 5 seconds after its onset.
                if elapsed < 0.5 + 5.0:
                    screen.blit(target_surface, target_rect)
            pygame.display.flip()
            
            # If trial duration (0.5 sec delay + 5 sec target) has elapsed, end trial.
            if elapsed >= 0.5 + 5.0:
                trial_done = True
            clock.tick(60)
        
        # ------------------------------------
        # Determine and record the response.
        # ------------------------------------
        if response is None:
            key_str = None
            # If no response is made, check if the correct answer was "none".
            if str(trial["corrAns"]).lower() == "none":
                corr = 1
            else:
                corr = 0
        else:
            # Map Pygame key to string.
            if response == pygame.K_LEFT:
                key_str = "left"
            elif response == pygame.K_RIGHT:
                key_str = "right"
            elif response == pygame.K_DOWN:
                key_str = "down"
            else:
                key_str = ""
            # Compare the response to the correct answer.
            corr = 1 if key_str == str(trial["corrAns"]) else 0
        
        # Append trial data.
        data_list.append({
            "trial": trial_number,
            "word": trial["word"],
            "colour": trial["colour"],
            "corrAns": trial["corrAns"],
            "response": key_str,
            "correct": corr,
            "RT": rt if rt is not None else ""
        })
        
        # Optional: brief inter-trial interval (e.g., 500 ms).
        pygame.time.delay(500)
    
    # =====================================
    # Save data to CSV file
    # =====================================
    try:
        with open(data_filename, "w", newline="") as csvfile:
            fieldnames = ["trial", "word", "colour", "corrAns", "response", "correct", "RT"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in data_list:
                writer.writerow(row)
        print("Data saved to", data_filename)
    except Exception as e:
        print("Error saving data:", e)
    
    # =====================================
    # End of Experiment Message
    # =====================================
    end_text = "Experiment finished. Press any key to exit."
    screen.fill((0, 0, 0))
    end_surface = instruction_font.render(end_text, True, (255, 255, 255))
    end_rect = end_surface.get_rect(center=(screen_width / 2, screen_height / 2))
    screen.blit(end_surface, end_rect)
    pygame.display.flip()
    
    waiting = True
    while waiting:
        for event in pygame.event.get():
            if event.type == pygame.KEYDOWN or event.type == pygame.QUIT:
                waiting = False
        clock.tick(60)
    
    pygame.quit()
    sys.exit()

if __name__ == "__main__":
    main()
