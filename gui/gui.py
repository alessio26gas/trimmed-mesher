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
        ctk.set_default_color_theme("theme.json")

        super().__init__(*args, **kwargs)
        self.title("Trimmed Mesher")
        self.resizable(True, True)
        self.geometry("540x300")
        self.minsize(540, 400)
        self.maxsize(640, 800)
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
        self.main_frame = ctk.CTkFrame(self, height=240, border_width=1, corner_radius=0)
        self.main_frame.grid_columnconfigure(0, weight=1)

        # Tab Buttons
        self.tabs_frame = ctk.CTkFrame(self.main_frame, height=36, border_width=1, corner_radius=0)
        self.tabs_frame.grid_columnconfigure([0, 1, 2, 3], weight=1)

        self.on = ['gray99', 'gray12']
        self.on_h = ["gray95","gray26"]
        self.off = ["gray90","gray20"]
        self.off_h = ["gray85","gray30"]

        self.general_button = ctk.CTkButton(
            self.tabs_frame, text="General Settings",
            command=lambda: self.switch_tab(0),
            corner_radius=0,
            border_width=0,
            fg_color=self.on, hover_color=self.on_h
        )
        self.general_button.grid(row=0, column=0, padx=(1, 0), pady=(1, 0), sticky="news")

        self.nwl_button = ctk.CTkButton(
            self.tabs_frame, text="Near wall cells",
            command=lambda: self.switch_tab(1),
            corner_radius=0,
            border_width=0,
            fg_color=self.off, hover_color=self.off_h
        )
        self.nwl_button.grid(row=0, column=1, padx=0, pady=(1, 0), sticky="news")

        self.smoothing_button = ctk.CTkButton(
            self.tabs_frame, text="Smoothing",
            command=lambda: self.switch_tab(2),
            corner_radius=0,
            border_width=0,
            fg_color=self.off, hover_color=self.off_h
        )
        self.smoothing_button.grid(row=0, column=2, padx=0, pady=(1, 0), sticky="news")

        self.coarsening_button = ctk.CTkButton(
            self.tabs_frame, text="Coarsening",
            command=lambda: self.switch_tab(3),
            corner_radius=0,
            border_width=0,
            fg_color=self.off, hover_color=self.off_h
        )
        self.coarsening_button.grid(row=0, column=3, padx=(0, 1), pady=(1, 0), sticky="news")

        self.tabs_frame.grid(row=0, column=0, padx=4, pady=(4, 0), sticky="news")

        # General Frame
        self.general_frame = ctk.CTkFrame(self.main_frame, height=200, border_width=1, corner_radius=0)

        ctk.CTkLabel(self.general_frame, text="Output file path").place(anchor="w", relx=0.01, rely=1/10)
        self.output_file_entry = ctk.CTkEntry(self.general_frame, corner_radius=0)
        self.output_file_entry.place(anchor="center", relwidth=0.5, relx=0.5, rely=1/10)
        ctk.CTkButton(
            self.general_frame, text="Browse...", corner_radius=0, command=self.output_file
        ).place(anchor="e", relwidth=0.24, relx=0.99, rely=1/10)

        ctk.CTkLabel(self.general_frame, text="Body curve path").place(anchor="w", relx=0.01, rely=3/10)
        self.curve_entry = ctk.CTkEntry(self.general_frame, corner_radius=0)
        self.curve_entry.place(anchor="center", relwidth=0.5, relx=0.5, rely=3/10)
        ctk.CTkButton(
            self.general_frame, text="Browse...", corner_radius=0, command=self.load_curve
        ).place(anchor="e", relwidth=0.24, relx=0.99, rely=3/10)

        ctk.CTkLabel(self.general_frame, text="Cell size").place(anchor="w", relx=0.01, rely=5/10)
        self.cell_size_entry = ctk.CTkEntry(self.general_frame)
        self.cell_size_entry.place(anchor="e", relwidth=0.24, relx=0.49, rely=5/10)

        ctk.CTkLabel(self.general_frame, text="Rows").place(anchor="w", relx=0.01, rely=7/10)
        self.rows_entry = ctk.CTkEntry(self.general_frame)
        self.rows_entry.place(anchor="e", relwidth=0.24, relx=0.49, rely=7/10)

        ctk.CTkLabel(self.general_frame, text="Columns").place(anchor="w", relx=0.01, rely=9/10)
        self.cols_entry = ctk.CTkEntry(self.general_frame)
        self.cols_entry.place(anchor="e", relwidth=0.24, relx=0.49, rely=9/10)

        ctk.CTkLabel(self.general_frame, text="Center (X,Y)").place(anchor="w", relx=0.51, rely=5/10)
        self.center_x_entry = ctk.CTkEntry(self.general_frame)
        self.center_x_entry.place(anchor="e", relwidth=0.12, relx=0.87, rely=5/10)
        self.center_y_entry = ctk.CTkEntry(self.general_frame)
        self.center_y_entry.place(anchor="e", relwidth=0.12, relx=0.99, rely=5/10)

        ctk.CTkLabel(self.general_frame, text="Rotation angle").place(anchor="w", relx=0.51, rely=7/10)
        self.rotation_angle_entry = ctk.CTkEntry(self.general_frame)
        self.rotation_angle_entry.place(anchor="e", relwidth=0.24, relx=0.99, rely=7/10)

        ctk.CTkLabel(self.general_frame, text="Rotation center (X,Y)").place(anchor="w", relx=0.51, rely=9/10)
        self.rotation_center_x_entry = ctk.CTkEntry(self.general_frame)
        self.rotation_center_x_entry.place(anchor="e", relwidth=0.12, relx=0.87, rely=9/10)
        self.rotation_center_y_entry = ctk.CTkEntry(self.general_frame)
        self.rotation_center_y_entry.place(anchor="e", relwidth=0.12, relx=0.99, rely=9/10)

        self.general_frame.grid(row=1, column=0, padx=4, pady=(0, 6), sticky="news")

        # Near Wall Cells Frame
        self.nwl_frame = ctk.CTkFrame(self.main_frame, height=200, border_width=1, corner_radius=0)

        ctk.CTkLabel(self.nwl_frame, text="Enable near wall layer").place(anchor="w", relx=0.01, rely=1/10)
        self.enable_nwl_cb = ctk.CTkCheckBox(self.nwl_frame, text="")
        self.enable_nwl_cb.place(anchor="w", relx=0.26, rely=1/10)

        ctk.CTkLabel(self.nwl_frame, text="Near wall thickness").place(anchor="w", relx=0.51, rely=1/10)
        self.nwl_first_entry = ctk.CTkEntry(self.nwl_frame)
        self.nwl_first_entry.place(anchor="e", relwidth=0.24, relx=0.99, rely=1/10)

        ctk.CTkLabel(self.nwl_frame, text="Total thickness").place(anchor="w", relx=0.01, rely=3/10)
        self.nwl_distance_entry = ctk.CTkEntry(self.nwl_frame)
        self.nwl_distance_entry.place(anchor="e", relwidth=0.24, relx=0.49, rely=3/10)

        ctk.CTkLabel(self.nwl_frame, text="Last cell thickness").place(anchor="w", relx=0.51, rely=3/10)
        self.nwl_last_entry = ctk.CTkEntry(self.nwl_frame)
        self.nwl_last_entry.place(anchor="e", relwidth=0.24, relx=0.99, rely=3/10)

        ctk.CTkLabel(self.nwl_frame, text="Number of Layers").place(anchor="w", relx=0.01, rely=5/10)
        self.nwl_n_entry = ctk.CTkEntry(self.nwl_frame)
        self.nwl_n_entry.place(anchor="e", relwidth=0.24, relx=0.49, rely=5/10)

        ctk.CTkLabel(self.nwl_frame, text="Stretch Factor").place(anchor="w", relx=0.51, rely=5/10)
        self.nwl_sf_entry = ctk.CTkEntry(self.nwl_frame)
        self.nwl_sf_entry.place(anchor="e", relwidth=0.24, relx=0.99, rely=5/10)

        self.nwl_frame.grid(row=1, column=0, padx=4, pady=(0, 6), sticky="news")

        # Smoothing Frame
        self.smoothing_frame = ctk.CTkFrame(self.main_frame, height=200, border_width=1, corner_radius=0)

        self.smoothing_frame.grid(row=1, column=0, padx=4, pady=(0, 6), sticky="news")

        # Coarsening Frame
        self.coarsening_frame = ctk.CTkFrame(self.main_frame, height=200, border_width=1, corner_radius=0)

        self.coarsening_frame.grid(row=1, column=0, padx=4, pady=(0, 6), sticky="news")

        # Main Frame Configuration
        self.general_frame.tkraise()
        self.main_frame.grid(row=0, column=0, columnspan=2, padx=4, pady=4, sticky="ew")

        # Control Center
        self.control_center = ctk.CTkFrame(self, height=32, border_width=1, corner_radius=0)
        self.control_center.grid_columnconfigure([0, 2], weight=1)
        self.control_center.grid_columnconfigure(1, weight=0)

        if self._get_appearance_mode() == "light":
            start_image = ctk.CTkImage(Image.open(resource_path("images/start.png")), size=(32, 32))
        else:
            start_image = ctk.CTkImage(Image.open(resource_path("images/start_dark.png")), size=(32, 32))

        self.start_button = ctk.CTkButton(
            self.control_center, 
            command=self.start,
            width=25,
            corner_radius=5,
            border_width=0,
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
            border_width=0,
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
            border_width=0,
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
            border_width=0,
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
        self.output_text.grid(row=0, column=0, padx=0, pady=0, sticky="nsew")

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
        self.input.__init__()
        self.input.curve = self.curve_entry.get()
        self.input.outputfile = self.output_file_entry.get()
        self.input.cell_size = self.cell_size_entry.get()
        self.input.rows = self.rows_entry.get()
        self.input.cols = self.cols_entry.get()
        self.input.center_x = self.center_x_entry.get()
        self.input.center_y = self.center_y_entry.get()
        self.input.rotation_angle = self.rotation_angle_entry.get()
        self.input.rotation_center_x = self.rotation_center_x_entry.get()
        self.input.rotation_center_y = self.rotation_center_y_entry.get()
        self.input.enable_nwl = self.enable_nwl_cb.get()
        self.input.nwl_first = self.nwl_first_entry.get()
        self.input.nwl_last = self.nwl_last_entry.get()
        self.input.nwl_distance = self.nwl_distance_entry.get()
        self.input.nwl_n = self.nwl_n_entry.get()
        self.input.nwl_SF = self.nwl_sf_entry.get()

    def load_curve(self):
        self.input.curve = filedialog.askopenfilename(title="Load curve file", filetypes=(("CSV files", "*.csv"), ("all files", "*.*")))
        self.curve_entry.delete("0", "end")
        self.curve_entry.insert("0", self.input.curve)
        self.curve_entry.xview_moveto(1.0)

    def output_file(self):
        self.input.outputfile = filedialog.asksaveasfilename(defaultextension=".msh")
        self.output_file_entry.delete("0", "end")
        self.output_file_entry.insert("0", self.input.outputfile)
        self.output_file_entry.xview_moveto(1.0)

    def switch_tab(self, tab):
        self.general_button.configure(fg_color=self.off, hover_color=self.off_h)
        self.nwl_button.configure(fg_color=self.off, hover_color=self.off_h)
        self.smoothing_button.configure(fg_color=self.off, hover_color=self.off_h)
        self.coarsening_button.configure(fg_color=self.off, hover_color=self.off_h)
        match tab:
            case 0:
                self.general_frame.tkraise()
                self.general_button.configure(fg_color=self.on, hover_color=self.on_h)
            case 1:
                self.nwl_frame.tkraise()
                self.nwl_button.configure(fg_color=self.on, hover_color=self.on_h)
            case 2:
                self.smoothing_frame.tkraise()
                self.smoothing_button.configure(fg_color=self.on, hover_color=self.on_h)
            case 3:
                self.coarsening_frame.tkraise()
                self.coarsening_button.configure(fg_color=self.on, hover_color=self.on_h)

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
                self.info, border_width=2, corner_radius=0, width=260, height=250,
            )
            self.infoframe.grid(column=0, row=0, padx=2, pady=2)

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
        self.cell_size = "0"
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
        self.rows = "0"
        self.cols = "0"
        self.center_x = "0.0"
        self.center_y = "0.0"
        self.rotation_angle = "0.0"
        self.rotation_center_x = "0.0"
        self.rotation_center_y = "0.0"
        self.smoothing = "0"
        self.smoothing_iterations = "0"
        self.enable_nwl = "0"
        self.nwl_first = "0.0"
        self.nwl_last = self.cell_size
        self.nwl_distance = "0.0"
        self.nwl_n = "0"
        self.nwl_SF = "0.0"
    
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