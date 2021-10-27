from tkinter import *

class Application(Frame):
    def createWidgets(self):
        #entryfield
        self.entry = Entry(self, width = 50)
        self.entry.pack({"side": "top"})
        self.entry.insert(0, "a b") #should be changed - just for now
        self.entry.focus()

        #button
        self.button = Button(self, width = 20, text = "Click")
        self.button.pack({"side": "top"})

        #canvas
        self.canvas = Canvas(self, bg = "yellow", width = 480, height = 440)

        #I put these here for now, but they will appear after buttonpress
        line1 = self.canvas.create_line(100, 100, 300, 100)
        line2 = self.canvas.create_line(100, 150, 300, 150)
        label1 = self.canvas.create_text(200, 80, text="a", fill="black", font=('Helvetica 15 bold'))
        label1 = self.canvas.create_text(200, 170, text="b", fill="black", font=('Helvetica 15 bold'))
        self.canvas.pack({"side": "bottom"})
        

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()
    
    #when clicking the button, the string is interpreted and drawn
    def fetch():
        pass
    
root = Tk()
root.title("RNADomainViz")
root.geometry('500x500')

app = Application(master=root)
app.mainloop()
root.destroy()