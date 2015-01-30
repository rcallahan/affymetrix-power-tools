////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

#pragma once

using namespace System;
using namespace System::ComponentModel;
using namespace System::Collections;
using namespace System::Windows::Forms;
using namespace System::Data;
using namespace System::Drawing;


namespace Midas
{
	/// <summary>
	/// Summary for MidasProcessingBox
	///
	/// WARNING: If you change the name of this class, you will need to change the
	///          'Resource File Name' property for the managed resource compiler tool
	///          associated with all .resx files this class depends on.  Otherwise,
	///          the designers will not be able to interact properly with localized
	///          resources associated with this form.
	/// </summary>
	public __gc class MidasProcessingBox : public System::Windows::Forms::Form
	{
	public:
		MidasProcessingBox(void)
		{
			InitializeComponent();
		}

	protected:
		void Dispose(Boolean disposing)
		{
			if (disposing && components)
			{
				components->Dispose();
			}
			__super::Dispose(disposing);
		}
	private: System::Windows::Forms::Label *  label1;
	private: System::Windows::Forms::PictureBox *  pictureBox1;

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container* components;

		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			System::Resources::ResourceManager *  resources = new System::Resources::ResourceManager(__typeof(Midas::MidasProcessingBox));
			this->label1 = new System::Windows::Forms::Label();
			this->pictureBox1 = new System::Windows::Forms::PictureBox();
			this->SuspendLayout();
			//
			// label1
			//
			this->label1->Location = System::Drawing::Point(24, 16);
			this->label1->Name = S"label1";
			this->label1->Size = System::Drawing::Size(112, 23);
			this->label1->TabIndex = 0;
			this->label1->Text = S"Performing analysis.";
			this->label1->TextAlign = System::Drawing::ContentAlignment::MiddleCenter;
			//
			// pictureBox1
			//
			this->pictureBox1->Image = (__try_cast<System::Drawing::Image *  >(resources->GetObject(S"pictureBox1.Image")));
			this->pictureBox1->Location = System::Drawing::Point(32, 56);
			this->pictureBox1->Name = S"pictureBox1";
			this->pictureBox1->Size = System::Drawing::Size(100, 32);
			this->pictureBox1->TabIndex = 1;
			this->pictureBox1->TabStop = false;
			//
			// MidasProcessingBox
			//
			this->AutoScaleBaseSize = System::Drawing::Size(5, 13);
			this->ClientSize = System::Drawing::Size(160, 102);
			this->ControlBox = false;
			this->Controls->Add(this->pictureBox1);
			this->Controls->Add(this->label1);
			this->MaximizeBox = false;
			this->MinimizeBox = false;
			this->Name = S"MidasProcessingBox";
			this->Text = S"Processing...";
			this->ResumeLayout(false);

		}
	};
}
