import customtkinter as ctk
import tkinter as tk
from tkinter import messagebox, filedialog
from process import ProcessHandler
from PIL import Image
from utils import get_font, resource_path

version = "dev"
copyright = "(C) 2025 Alessio Improta"


class TrimmedMesher(ctk.CTk):
    def __init__(self, *args, **kwargs):
        ctk.set_appearance_mode("system")

        super().__init__(*args, **kwargs)
        self.title("Trimmed Mesher")
        self.resizable(True, True)
        self.geometry("540x300")
        self.minsize(540, 350)
        self.maxsize(640, 720)
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(2, weight=1)
        self.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.after(
            201,
            lambda: self.iconphoto(
                False, tk.PhotoImage(file=resource_path("images/icon.png"))
            )
        )

        # Process Handler
        self.process_handler = ProcessHandler(self.append_output)
        self.input = Input()
        self.info = None

        # Main Frame
        self.main_frame = ctk.CTkFrame(self, height=200, border_width=1, corner_radius=0)

        self.curve_button = ctk.CTkButton(self.main_frame, text="Load curve", command=self.load_curve)
        self.curve_button.place(anchor="center", relx=0.5, rely=0.1)

        self.output_button = ctk.CTkButton(self.main_frame, text="Output file", command=self.output_file)
        self.output_button.place(anchor="center", relx=0.5, rely=0.3)

        self.input_entry = ctk.CTkEntry(self.main_frame, placeholder_text="Input field")
        self.input_entry.place(anchor="center", relx=0.5, rely=0.5)

        self.main_frame.grid(row=0, column=0, columnspan=2, padx=4, pady=4, sticky="ew")

        # Control Center
        self.control_center = ctk.CTkFrame(self, height=32, border_width=1, corner_radius=0)
        self.control_center.grid_columnconfigure(0, weight=1)
        self.control_center.grid_columnconfigure(1, weight=0)
        self.control_center.grid_columnconfigure(2, weight=1)

        if self._get_appearance_mode() == "light":
            start_image = ctk.CTkImage(Image.open(resource_path("images/start.png")), size=(32, 32))
        else:
            start_image = ctk.CTkImage(Image.open(resource_path("images/start_dark.png")), size=(32, 32))

        self.start_button = ctk.CTkButton(
            self.control_center, 
            command=self.start,
            width=25,
            corner_radius=5,
            text="",
            image=start_image,
            fg_color="transparent",
            text_color=["gray10", "gray95"],
            hover_color=["gray95", "gray10"]
        )
        self.start_button.grid(row=0, column=0, sticky="e", padx=(0, 10), pady=4)

        if self._get_appearance_mode() == "light":
            stop_image = ctk.CTkImage(Image.open(resource_path("images/stop.png")), size=(32, 32))
        else:
            stop_image = ctk.CTkImage(Image.open(resource_path("images/stop_dark.png")), size=(32, 32))

        self.stop_button = ctk.CTkButton(
            self.control_center, 
            command=self.stop,
            width=25,
            corner_radius=5,
            text="",
            image=stop_image,
            fg_color="transparent",
            text_color=["gray10", "gray95"],
            hover_color=["gray95", "gray10"]
        )
        self.stop_button.grid(row=0, column=2, sticky="w", padx=(10, 0), pady=4)

        if self._get_appearance_mode() == "light":
            appearance_image = ctk.CTkImage(Image.open(resource_path("images/moon.png")), size=(20, 20))
        else:
            appearance_image = ctk.CTkImage(Image.open(resource_path("images/sun.png")), size=(20, 20))

        self.appearance_button = ctk.CTkButton(
            self.control_center,
            command=self.toggle_appearance_mode,
            width=25,
            corner_radius=5,
            text="",
            image=appearance_image,
            fg_color="transparent",
            text_color=["gray10", "gray95"],
            hover_color=["gray95", "gray10"]
        )
        self.appearance_button.place(anchor="w", relx=0.01, rely=0.5)

        if self._get_appearance_mode() == "light":
            info_image = ctk.CTkImage(Image.open(resource_path("images/info.png")), size=(16, 16))
        else:
            info_image = ctk.CTkImage(Image.open(resource_path("images/info_dark.png")), size=(16, 16))

        self.info_button = ctk.CTkButton(
            self.control_center,
            command=self.get_info,
            width=25,
            corner_radius=5,
            text="",
            image=info_image,
            fg_color="transparent",
            text_color=["gray10", "gray95"],
            hover_color=["gray95", "gray10"]
        )
        self.info_button.place(anchor="e", relx=0.99, rely=0.5)

        self.control_center.grid(row=1, column=0, columnspan=2, padx=4, pady=0, sticky="ew")

        # Command Window
        self.command_window = ctk.CTkFrame(self, border_width=1, corner_radius=0)
        self.command_window.grid_rowconfigure(0, weight=1)
        self.command_window.grid_columnconfigure(0, weight=1)

        self.output_text = ctk.CTkTextbox(
            self.command_window,
            state="disabled",
            wrap="none",
            corner_radius=0,
            border_width=2,
            font=get_font(),
        )
        self.output_text.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        self.command_window.grid(row=2, column=0, columnspan=2, padx=4, pady=4, sticky="nsew")

    def start(self):
        self.get_input()
        self.process_handler.start_process(self.input.split())

    def stop(self):
        self.process_handler.stop_process()

    def append_output(self, text):
        self.output_text.configure(state="normal")
        self.output_text.insert("end", text)
        self.output_text.yview_moveto(0.99)
        self.output_text.configure(state="disabled")

    def get_input(self):
        self.input.cell_size = self.input_entry.get()

    def load_curve(self):
        self.input.curve = filedialog.askopenfilename(title="Load curve file", filetypes=(("CSV files", "*.csv"), ("all files", "*.*")))
    
    def output_file(self):
        self.input.outputfile = filedialog.asksaveasfilename(defaultextension=".msh")

    def get_info(self):
        if self.info is None or not self.info.winfo_exists():
            self.info = ctk.CTkToplevel()
            self.info.title("About")
            self.info.configure(width=260, height=250)
            self.info.resizable(False, False)
            self.info.after(
                201,
                lambda: self.info.iconphoto(
                    False, tk.PhotoImage(file=resource_path("images/icon.png"))
                ),
            )

            self.infoframe = ctk.CTkFrame(
                self.info, border_width=3, corner_radius=0, width=260, height=250,
            )
            self.infoframe.grid(column=0, row=0)

            icon = ctk.CTkImage(Image.open(resource_path("images/icon.png")), size=(128, 128))
            ctk.CTkLabel(
                self.infoframe, text="", image=icon
            ).place(anchor="center", relx=0.5, y=70)
            ctk.CTkLabel(
                self.infoframe, text="Trimmed Mesher", font=("Sans", 20)
            ).place(anchor="center", relx=0.5, y=140)
            ctk.CTkLabel(
                self.infoframe, text="Version " + version
            ).place(anchor="center", relx=0.5, y=163)
            ctk.CTkLabel(
                self.infoframe, text=copyright
            ).place(anchor="center", relx=0.5, y=185)
            ctk.CTkLabel(
                self.infoframe, text=(
                "Main author: Alessio Improta\n"
                + "alessio.improta@polito.it"
                )
            ).place(anchor="center", relx=0.5, rely=0.88)

            self.info.after(50, self.info.lift)
            self.info.after(50, self.info.focus)

        else:
            self.info.lift()
            self.info.focus()
    
    def toggle_appearance_mode(self):
        if self._get_appearance_mode() == "dark":
            new_mode = "light"
            self.appearance_button.configure(image=ctk.CTkImage(Image.open(resource_path("images/moon.png")), size=(20, 20)))
            self.info_button.configure(image=ctk.CTkImage(Image.open(resource_path("images/info.png")), size=(16, 16)))
            self.start_button.configure(image=ctk.CTkImage(Image.open(resource_path("images/start.png")), size=(32, 32)))
            self.stop_button.configure(image=ctk.CTkImage(Image.open(resource_path("images/stop.png")), size=(32, 32)))
        else:
            new_mode = "dark"
            self.appearance_button.configure(image=ctk.CTkImage(Image.open(resource_path("images/sun.png")), size=(20, 20)))
            self.info_button.configure(image=ctk.CTkImage(Image.open(resource_path("images/info_dark.png")), size=(16, 16)))
            self.start_button.configure(image=ctk.CTkImage(Image.open(resource_path("images/start_dark.png")), size=(32, 32)))
            self.stop_button.configure(image=ctk.CTkImage(Image.open(resource_path("images/stop_dark.png")), size=(32, 32)))
        ctk.set_appearance_mode(new_mode)

    def on_closing(self):
        quit_ = messagebox.askokcancel(title="Exit?", message="Do you want to exit?")
        if quit_:
            self.process_handler.cleanup()
            self.destroy()


class Input():
    def __init__(self):
        self.curve = "0"
        self.outputfile = "mesh.msh"
        self.cell_size = "0.01"
        self.coarsening_level_0 = "0"
        self.coarsening_level_1 = "0"
        self.coarsening_level_2 = "0"
        self.coarsening_level_3 = "0"
        self.coarsening_cells_0 = "0"
        self.coarsening_cells_1 = "0"
        self.coarsening_cells_2 = "0"
        self.coarsening_cells_3 = "0"
        self.fast_coarsening = "0"
        self.conformal_coarsening = "0"
        self.rows = "256"
        self.cols = "256"
        self.center_x = "0.0"
        self.center_y = "0.0"
        self.rotation_angle = "0.0"
        self.rotation_center_x = "0.0"
        self.rotation_center_y = "0.0"
        self.smoothing = "1"
        self.smoothing_iterations = "20000"
        self.enable_nwl = "1"
        self.nwl_first = "0.0001"
        self.nwl_last = self.cell_size
        self.nwl_distance = "0.0"
        self.nwl_n = "0"
        self.nwl_SF = "1.25"
    
    def split(self):
        return [
            self.curve,
            self.outputfile,
            self.cell_size,
            self.coarsening_level_0,
            self.coarsening_level_1,
            self.coarsening_level_2,
            self.coarsening_level_3,
            self.coarsening_cells_0,
            self.coarsening_cells_1,
            self.coarsening_cells_2,
            self.coarsening_cells_3,
            self.fast_coarsening,
            self.conformal_coarsening,
            self.rows,
            self.cols,
            self.center_x,
            self.center_y,
            self.rotation_angle,
            self.rotation_center_x,
            self.rotation_center_y,
            self.smoothing,
            self.smoothing_iterations,
            self.enable_nwl,
            self.nwl_first,
            self.nwl_last,
            self.nwl_distance,
            self.nwl_n,
            self.nwl_SF
        ]


if __name__ == "__main__":
    TrimmedMesher(className="Trimmed Mesher").mainloop()