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

#include <string>
#using <mscorlib.dll>
#include <stdexcept>
#include "MidasConfigureRun.h"
#include "FolderDialog.h"
// want to use MessageBox from System::Windows::Forms,
// not the SDK MessageBox()
#ifdef MessageBox
#undef MessageBox
#endif
// similarly for GetObject
#ifdef GetObject
#undef GetObject
#endif
#include "MidasProcessingBox.h"

namespace Midas
{
	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	using namespace System::Runtime::InteropServices;
	using namespace System::IO;

	/// <summary>
	/// Summary for WinMidas
	///
	/// WARNING: If you change the name of this class, you will need to change the
	///          'Resource File Name' property for the managed resource compiler tool
	///          associated with all .resx files this class depends on.  Otherwise,
	///          the designers will not be able to interact properly with localized
	///          resources associated with this form.
	/// </summary>
	public __gc class WinMidas : public System::Windows::Forms::Form
	{
	public:
		WinMidas(const std::string& execVersion)
		: execVersion(execVersion)
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
	private: const std::string& execVersion;
	private: System::Windows::Forms::OpenFileDialog *  openFileDialog1;
	private: System::Windows::Forms::MainMenu *  mainMenu1;


	private: System::Windows::Forms::Label *  label1;
	private: System::Windows::Forms::TextBox *  celsFileNameTextBox;
	private: System::Windows::Forms::Button *  celsFileButton;
	private: System::Windows::Forms::Label *  label2;

	private: System::Windows::Forms::Button *  geneDataFileButton;
	private: System::Windows::Forms::TextBox *  geneDataFileTextBox;
	private: System::Windows::Forms::Label *  label3;
	private: System::Windows::Forms::TextBox *  probesetDataFileTextBox;
	private: System::Windows::Forms::Button *  probesetDataFileButton;
	private: System::Windows::Forms::Label *  label4;
	private: System::Windows::Forms::TextBox *  mapFileTextBox;
	private: System::Windows::Forms::Button *  mapFileButton;
	private: System::Windows::Forms::Label *  label5;
	private: System::Windows::Forms::TextBox *  outFolderTextBox;


	private: System::Windows::Forms::Label *  label6;
	private: System::Windows::Forms::TextBox *  logStabilizeTextBox;
	private: System::Windows::Forms::GroupBox *  pValuesGroupBox;
	private: System::Windows::Forms::RadioButton *  pvaluesRadioYes;
	private: System::Windows::Forms::RadioButton *  pvaluesRadioNo;

	private: bool notWantPvalues;
	private: bool wantFstats;
	private: bool wantIndices;
	private: bool noLogTransform;
	private: System::Windows::Forms::GroupBox *  fstatsGroupBox;
	private: System::Windows::Forms::RadioButton *  fstatsRadioYes;
	private: System::Windows::Forms::RadioButton *  fstatsRadioNo;
	private: System::Windows::Forms::GroupBox *  indicesGroupBox;
	private: System::Windows::Forms::RadioButton *  indicesRadioYes;
	private: System::Windows::Forms::RadioButton *  indicesRadioNo;
	private: System::Windows::Forms::Button *  analyzeButton;
	private: System::Windows::Forms::Button *  exitButton;
	private: System::Windows::Forms::Label *  label7;
	private: System::Windows::Forms::Button *  outFolderButton;
	private: System::Windows::Forms::GroupBox *  logTransGroupBox;
	private: System::Windows::Forms::RadioButton *  logTransRadioYes;
	private: System::Windows::Forms::RadioButton *  logTransRadioNo;




	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container * components;

		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->openFileDialog1 = new System::Windows::Forms::OpenFileDialog();
			this->mainMenu1 = new System::Windows::Forms::MainMenu();
			this->label1 = new System::Windows::Forms::Label();
			this->celsFileNameTextBox = new System::Windows::Forms::TextBox();
			this->celsFileButton = new System::Windows::Forms::Button();
			this->label2 = new System::Windows::Forms::Label();
			this->geneDataFileTextBox = new System::Windows::Forms::TextBox();
			this->geneDataFileButton = new System::Windows::Forms::Button();
			this->label3 = new System::Windows::Forms::Label();
			this->probesetDataFileTextBox = new System::Windows::Forms::TextBox();
			this->probesetDataFileButton = new System::Windows::Forms::Button();
			this->label4 = new System::Windows::Forms::Label();
			this->mapFileTextBox = new System::Windows::Forms::TextBox();
			this->mapFileButton = new System::Windows::Forms::Button();
			this->label5 = new System::Windows::Forms::Label();
			this->outFolderTextBox = new System::Windows::Forms::TextBox();
			this->label6 = new System::Windows::Forms::Label();
			this->logStabilizeTextBox = new System::Windows::Forms::TextBox();
			this->pValuesGroupBox = new System::Windows::Forms::GroupBox();
			this->pvaluesRadioNo = new System::Windows::Forms::RadioButton();
			this->pvaluesRadioYes = new System::Windows::Forms::RadioButton();
			this->fstatsGroupBox = new System::Windows::Forms::GroupBox();
			this->fstatsRadioNo = new System::Windows::Forms::RadioButton();
			this->fstatsRadioYes = new System::Windows::Forms::RadioButton();
			this->indicesGroupBox = new System::Windows::Forms::GroupBox();
			this->indicesRadioNo = new System::Windows::Forms::RadioButton();
			this->indicesRadioYes = new System::Windows::Forms::RadioButton();
			this->analyzeButton = new System::Windows::Forms::Button();
			this->exitButton = new System::Windows::Forms::Button();
			this->label7 = new System::Windows::Forms::Label();
			this->outFolderButton = new System::Windows::Forms::Button();
			this->logTransGroupBox = new System::Windows::Forms::GroupBox();
			this->logTransRadioNo = new System::Windows::Forms::RadioButton();
			this->logTransRadioYes = new System::Windows::Forms::RadioButton();
			this->pValuesGroupBox->SuspendLayout();
			this->fstatsGroupBox->SuspendLayout();
			this->indicesGroupBox->SuspendLayout();
			this->logTransGroupBox->SuspendLayout();
			this->SuspendLayout();
			//
			// label1
			//
			this->label1->Location = System::Drawing::Point(16, 64);
			this->label1->Name = S"label1";
			this->label1->Size = System::Drawing::Size(152, 23);
			this->label1->TabIndex = 0;
			this->label1->Text = S"Cel Group File:";
			this->label1->TextAlign = System::Drawing::ContentAlignment::MiddleRight;
			this->label1->Click += new System::EventHandler(this, &Midas::WinMidas::label1_Click);
			//
			// celsFileNameTextBox
			//
			this->celsFileNameTextBox->Location = System::Drawing::Point(176, 64);
			this->celsFileNameTextBox->Name = S"celsFileNameTextBox";
			this->celsFileNameTextBox->Size = System::Drawing::Size(264, 20);
			this->celsFileNameTextBox->TabIndex = 1;
			this->celsFileNameTextBox->Text = S"";
			//
			// celsFileButton
			//
			this->celsFileButton->Location = System::Drawing::Point(456, 64);
			this->celsFileButton->Name = S"celsFileButton";
			this->celsFileButton->TabIndex = 2;
			this->celsFileButton->Text = S"Select...";
			this->celsFileButton->Click += new System::EventHandler(this, &Midas::WinMidas::celsFileButton_Click);
			//
			// label2
			//
			this->label2->Location = System::Drawing::Point(16, 104);
			this->label2->Name = S"label2";
			this->label2->Size = System::Drawing::Size(152, 23);
			this->label2->TabIndex = 3;
			this->label2->Text = S"Gene Level Summary File:";
			this->label2->TextAlign = System::Drawing::ContentAlignment::MiddleRight;
			//
			// geneDataFileTextBox
			//
			this->geneDataFileTextBox->Location = System::Drawing::Point(176, 104);
			this->geneDataFileTextBox->Name = S"geneDataFileTextBox";
			this->geneDataFileTextBox->Size = System::Drawing::Size(264, 20);
			this->geneDataFileTextBox->TabIndex = 4;
			this->geneDataFileTextBox->Text = S"";
			//
			// geneDataFileButton
			//
			this->geneDataFileButton->Location = System::Drawing::Point(456, 104);
			this->geneDataFileButton->Name = S"geneDataFileButton";
			this->geneDataFileButton->TabIndex = 5;
			this->geneDataFileButton->Text = S"Select...";
			this->geneDataFileButton->Click += new System::EventHandler(this, &Midas::WinMidas::geneDataFileButton_Click);
			//
			// label3
			//
			this->label3->Location = System::Drawing::Point(16, 144);
			this->label3->Name = S"label3";
			this->label3->Size = System::Drawing::Size(152, 23);
			this->label3->TabIndex = 6;
			this->label3->Text = S"Exon Level Summary File:";
			this->label3->TextAlign = System::Drawing::ContentAlignment::MiddleRight;
			//
			// probesetDataFileTextBox
			//
			this->probesetDataFileTextBox->Location = System::Drawing::Point(176, 144);
			this->probesetDataFileTextBox->Name = S"probesetDataFileTextBox";
			this->probesetDataFileTextBox->Size = System::Drawing::Size(264, 20);
			this->probesetDataFileTextBox->TabIndex = 7;
			this->probesetDataFileTextBox->Text = S"";
			//
			// probesetDataFileButton
			//
			this->probesetDataFileButton->Location = System::Drawing::Point(456, 144);
			this->probesetDataFileButton->Name = S"probesetDataFileButton";
			this->probesetDataFileButton->TabIndex = 8;
			this->probesetDataFileButton->Text = S"Select...";
			this->probesetDataFileButton->Click += new System::EventHandler(this, &Midas::WinMidas::probesetDataFileButton_Click);
			//
			// label4
			//
			this->label4->Location = System::Drawing::Point(16, 184);
			this->label4->Name = S"label4";
			this->label4->Size = System::Drawing::Size(152, 23);
			this->label4->TabIndex = 9;
			this->label4->Text = S"Meta Probeset List:";
			this->label4->TextAlign = System::Drawing::ContentAlignment::MiddleRight;
			//
			// mapFileTextBox
			//
			this->mapFileTextBox->Location = System::Drawing::Point(176, 184);
			this->mapFileTextBox->Name = S"mapFileTextBox";
			this->mapFileTextBox->Size = System::Drawing::Size(264, 20);
			this->mapFileTextBox->TabIndex = 10;
			this->mapFileTextBox->Text = S"";
			//
			// mapFileButton
			//
			this->mapFileButton->Location = System::Drawing::Point(456, 184);
			this->mapFileButton->Name = S"mapFileButton";
			this->mapFileButton->TabIndex = 11;
			this->mapFileButton->Text = S"Select...";
			this->mapFileButton->Click += new System::EventHandler(this, &Midas::WinMidas::mapFileButton_Click);
			//
			// label5
			//
			this->label5->Location = System::Drawing::Point(16, 224);
			this->label5->Name = S"label5";
			this->label5->Size = System::Drawing::Size(152, 23);
			this->label5->TabIndex = 12;
			this->label5->Text = S"Output Folder:";
			this->label5->TextAlign = System::Drawing::ContentAlignment::MiddleRight;
			//
			// outFolderTextBox
			//
			this->outFolderTextBox->Location = System::Drawing::Point(176, 224);
			this->outFolderTextBox->Name = S"outFolderTextBox";
			this->outFolderTextBox->Size = System::Drawing::Size(264, 20);
			this->outFolderTextBox->TabIndex = 13;
			this->outFolderTextBox->Text = S"";
			//
			// label6
			//
			this->label6->Location = System::Drawing::Point(200, 264);
			this->label6->Name = S"label6";
			this->label6->Size = System::Drawing::Size(128, 23);
			this->label6->TabIndex = 15;
			this->label6->Text = S"Log Stabilization Factor:";
			this->label6->TextAlign = System::Drawing::ContentAlignment::MiddleRight;
			//
			// logStabilizeTextBox
			//
			this->logStabilizeTextBox->Location = System::Drawing::Point(352, 264);
			this->logStabilizeTextBox->Name = S"logStabilizeTextBox";
			this->logStabilizeTextBox->Size = System::Drawing::Size(48, 20);
			this->logStabilizeTextBox->TabIndex = 16;
			this->logStabilizeTextBox->Text = S"8.0";
			this->logStabilizeTextBox->TextAlign = System::Windows::Forms::HorizontalAlignment::Center;
			//
			// pValuesGroupBox
			//
			this->pValuesGroupBox->Controls->Add(this->pvaluesRadioNo);
			this->pValuesGroupBox->Controls->Add(this->pvaluesRadioYes);
			this->pValuesGroupBox->Location = System::Drawing::Point(48, 296);
			this->pValuesGroupBox->Name = S"pValuesGroupBox";
			this->pValuesGroupBox->Size = System::Drawing::Size(104, 72);
			this->pValuesGroupBox->TabIndex = 17;
			this->pValuesGroupBox->TabStop = false;
			this->pValuesGroupBox->Text = S"Save P Values";
			this->pValuesGroupBox->Enter += new System::EventHandler(this, &Midas::WinMidas::pValuesGroupBox_Enter);
			//
			// pvaluesRadioNo
			//
			this->pvaluesRadioNo->Location = System::Drawing::Point(24, 40);
			this->pvaluesRadioNo->Name = S"pvaluesRadioNo";
			this->pvaluesRadioNo->Size = System::Drawing::Size(64, 24);
			this->pvaluesRadioNo->TabIndex = 1;
			this->pvaluesRadioNo->Text = S"No";
			this->pvaluesRadioNo->Click += new System::EventHandler(this, &Midas::WinMidas::pvaluesRadioYes_Click);
			//
			// pvaluesRadioYes
			//
			this->pvaluesRadioYes->Checked = true;
			this->pvaluesRadioYes->Location = System::Drawing::Point(24, 16);
			this->pvaluesRadioYes->Name = S"pvaluesRadioYes";
			this->pvaluesRadioYes->Size = System::Drawing::Size(56, 24);
			this->pvaluesRadioYes->TabIndex = 0;
			this->pvaluesRadioYes->TabStop = true;
			this->pvaluesRadioYes->Text = S"Yes";
			this->pvaluesRadioYes->Click += new System::EventHandler(this, &Midas::WinMidas::pvaluesRadioYes_Click);
			//
			// fstatsGroupBox
			//
			this->fstatsGroupBox->Controls->Add(this->fstatsRadioNo);
			this->fstatsGroupBox->Controls->Add(this->fstatsRadioYes);
			this->fstatsGroupBox->Location = System::Drawing::Point(168, 296);
			this->fstatsGroupBox->Name = S"fstatsGroupBox";
			this->fstatsGroupBox->Size = System::Drawing::Size(112, 72);
			this->fstatsGroupBox->TabIndex = 18;
			this->fstatsGroupBox->TabStop = false;
			this->fstatsGroupBox->Text = S"Save F Statistics";
			//
			// fstatsRadioNo
			//
			this->fstatsRadioNo->Checked = true;
			this->fstatsRadioNo->Location = System::Drawing::Point(32, 40);
			this->fstatsRadioNo->Name = S"fstatsRadioNo";
			this->fstatsRadioNo->Size = System::Drawing::Size(56, 24);
			this->fstatsRadioNo->TabIndex = 1;
			this->fstatsRadioNo->TabStop = true;
			this->fstatsRadioNo->Text = S"No";
			this->fstatsRadioNo->Click += new System::EventHandler(this, &Midas::WinMidas::fstatsRadioYes_Click);
			//
			// fstatsRadioYes
			//
			this->fstatsRadioYes->Location = System::Drawing::Point(32, 16);
			this->fstatsRadioYes->Name = S"fstatsRadioYes";
			this->fstatsRadioYes->Size = System::Drawing::Size(48, 24);
			this->fstatsRadioYes->TabIndex = 0;
			this->fstatsRadioYes->Text = S"Yes";
			this->fstatsRadioYes->Click += new System::EventHandler(this, &Midas::WinMidas::fstatsRadioYes_Click);
			//
			// indicesGroupBox
			//
			this->indicesGroupBox->Controls->Add(this->indicesRadioNo);
			this->indicesGroupBox->Controls->Add(this->indicesRadioYes);
			this->indicesGroupBox->Location = System::Drawing::Point(296, 296);
			this->indicesGroupBox->Name = S"indicesGroupBox";
			this->indicesGroupBox->Size = System::Drawing::Size(144, 72);
			this->indicesGroupBox->TabIndex = 19;
			this->indicesGroupBox->TabStop = false;
			this->indicesGroupBox->Text = S"Save Normalized Signal";
			this->indicesGroupBox->Enter += new System::EventHandler(this, &Midas::WinMidas::indicesGroupBox_Enter);
			//
			// indicesRadioNo
			//
			this->indicesRadioNo->Checked = true;
			this->indicesRadioNo->Location = System::Drawing::Point(40, 40);
			this->indicesRadioNo->Name = S"indicesRadioNo";
			this->indicesRadioNo->Size = System::Drawing::Size(48, 24);
			this->indicesRadioNo->TabIndex = 1;
			this->indicesRadioNo->TabStop = true;
			this->indicesRadioNo->Text = S"No";
			this->indicesRadioNo->Click += new System::EventHandler(this, &Midas::WinMidas::indicesRadioYes_Click);
			//
			// indicesRadioYes
			//
			this->indicesRadioYes->Location = System::Drawing::Point(40, 16);
			this->indicesRadioYes->Name = S"indicesRadioYes";
			this->indicesRadioYes->Size = System::Drawing::Size(48, 24);
			this->indicesRadioYes->TabIndex = 0;
			this->indicesRadioYes->Text = S"Yes";
			this->indicesRadioYes->Click += new System::EventHandler(this, &Midas::WinMidas::indicesRadioYes_Click);
			//
			// analyzeButton
			//
			this->analyzeButton->Location = System::Drawing::Point(168, 392);
			this->analyzeButton->Name = S"analyzeButton";
			this->analyzeButton->Size = System::Drawing::Size(96, 23);
			this->analyzeButton->TabIndex = 20;
			this->analyzeButton->Text = S"Analyze";
			this->analyzeButton->Click += new System::EventHandler(this, &Midas::WinMidas::analyzeButton_Click);
			//
			// exitButton
			//
			this->exitButton->Location = System::Drawing::Point(368, 392);
			this->exitButton->Name = S"exitButton";
			this->exitButton->Size = System::Drawing::Size(88, 23);
			this->exitButton->TabIndex = 21;
			this->exitButton->Text = S"Exit";
			this->exitButton->Click += new System::EventHandler(this, &Midas::WinMidas::exitButton_Click);
			//
			// label7
			//
			this->label7->Location = System::Drawing::Point(136, 8);
			this->label7->Name = S"label7";
			this->label7->Size = System::Drawing::Size(368, 48);
			this->label7->TabIndex = 22;
			this->label7->Text = S"Please select a cel group file, gene and exon level summary files, a meta probese"
				S"t list, and an output folder, and click on analyze.";
			this->label7->TextAlign = System::Drawing::ContentAlignment::MiddleCenter;
			//
			// outFolderButton
			//
			this->outFolderButton->Location = System::Drawing::Point(456, 224);
			this->outFolderButton->Name = S"outFolderButton";
			this->outFolderButton->TabIndex = 23;
			this->outFolderButton->Text = S"Select...";
			this->outFolderButton->Click += new System::EventHandler(this, &Midas::WinMidas::outFolderButton_Click);
			//
			// logTransGroupBox
			//
			this->logTransGroupBox->Controls->Add(this->logTransRadioNo);
			this->logTransGroupBox->Controls->Add(this->logTransRadioYes);
			this->logTransGroupBox->Location = System::Drawing::Point(456, 296);
			this->logTransGroupBox->Name = S"logTransGroupBox";
			this->logTransGroupBox->Size = System::Drawing::Size(128, 72);
			this->logTransGroupBox->TabIndex = 24;
			this->logTransGroupBox->TabStop = false;
			this->logTransGroupBox->Text = S"Log Transform Data";
			//
			// logTransRadioNo
			//
			this->logTransRadioNo->Location = System::Drawing::Point(32, 40);
			this->logTransRadioNo->Name = S"logTransRadioNo";
			this->logTransRadioNo->Size = System::Drawing::Size(56, 24);
			this->logTransRadioNo->TabIndex = 1;
			this->logTransRadioNo->Text = S"No";
			this->logTransRadioNo->Click += new System::EventHandler(this, &Midas::WinMidas::logTransRadioNo_Click);
			//
			// logTransRadioYes
			//
			this->logTransRadioYes->Checked = true;
			this->logTransRadioYes->Location = System::Drawing::Point(32, 16);
			this->logTransRadioYes->Name = S"logTransRadioYes";
			this->logTransRadioYes->Size = System::Drawing::Size(56, 24);
			this->logTransRadioYes->TabIndex = 0;
			this->logTransRadioYes->TabStop = true;
			this->logTransRadioYes->Text = S"Yes";
			this->logTransRadioYes->Click += new System::EventHandler(this, &Midas::WinMidas::logTransRadioNo_Click);
			//
			// WinMidas
			//
			this->AutoScaleBaseSize = System::Drawing::Size(5, 13);
			this->ClientSize = System::Drawing::Size(616, 441);
			this->Controls->Add(this->logTransGroupBox);
			this->Controls->Add(this->outFolderButton);
			this->Controls->Add(this->label7);
			this->Controls->Add(this->exitButton);
			this->Controls->Add(this->analyzeButton);
			this->Controls->Add(this->indicesGroupBox);
			this->Controls->Add(this->fstatsGroupBox);
			this->Controls->Add(this->pValuesGroupBox);
			this->Controls->Add(this->logStabilizeTextBox);
			this->Controls->Add(this->outFolderTextBox);
			this->Controls->Add(this->mapFileTextBox);
			this->Controls->Add(this->probesetDataFileTextBox);
			this->Controls->Add(this->geneDataFileTextBox);
			this->Controls->Add(this->celsFileNameTextBox);
			this->Controls->Add(this->label6);
			this->Controls->Add(this->label5);
			this->Controls->Add(this->mapFileButton);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->probesetDataFileButton);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->geneDataFileButton);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->celsFileButton);
			this->Controls->Add(this->label1);
			this->Menu = this->mainMenu1;
			this->Name = S"WinMidas";
			this->Text = S"MiDAS Alternative Splicing Analysis";
			this->pValuesGroupBox->ResumeLayout(false);
			this->fstatsGroupBox->ResumeLayout(false);
			this->indicesGroupBox->ResumeLayout(false);
			this->logTransGroupBox->ResumeLayout(false);
			this->ResumeLayout(false);

		}
	private: System::Void menuItem2_Click(System::Object *  sender, System::EventArgs *  e)
			 {
				 Application::Exit();
			 }

	private: System::Void celsFileButton_Click(System::Object *  sender, System::EventArgs *  e)
			 {
				 if(openFileDialog1->ShowDialog() == DialogResult::OK)
					 celsFileNameTextBox->Text = openFileDialog1->FileName;
			 }

private: System::Void geneDataFileButton_Click(System::Object *  sender, System::EventArgs *  e)
		 {
			 if (openFileDialog1->ShowDialog() == DialogResult::OK)
				 geneDataFileTextBox->Text = openFileDialog1->FileName;
		 }

private: System::Void probesetDataFileButton_Click(System::Object *  sender, System::EventArgs *  e)
		 {
			 if (openFileDialog1->ShowDialog() == DialogResult::OK)
				 probesetDataFileTextBox->Text = openFileDialog1->FileName;
		 }

private: System::Void mapFileButton_Click(System::Object *  sender, System::EventArgs *  e)
		 {
			 if (openFileDialog1->ShowDialog() == DialogResult::OK)
				 mapFileTextBox->Text = openFileDialog1->FileName;
		 }


private: System::Void pvaluesRadioYes_Click(System::Object *  sender, System::EventArgs *  e)
		 {
			 // use non-intuitive notWantPvalues to avoid problems
			 // with the initialization of bools in InitializeComponent()
			 if (sender == pvaluesRadioYes)
				 this->notWantPvalues = false;
			 else if (sender == pvaluesRadioNo)
				 this->notWantPvalues = true;
		 }

private: System::Void fstatsRadioYes_Click(System::Object *  sender, System::EventArgs *  e)
		 {
			 if (sender == fstatsRadioYes)
				 this->wantFstats = true;
			 else if (sender == fstatsRadioNo)
				 this->wantFstats = false;
		 }

private: System::Void indicesRadioYes_Click(System::Object *  sender, System::EventArgs *  e)
		 {
			 if (sender == indicesRadioYes)
				 this->wantIndices = true;
			 else if (sender == indicesRadioNo)
				 this->wantIndices = false;
		 }

private: System::Void analyzeButton_Click(System::Object *  sender, System::EventArgs *  e)
		 {
			 // require cels, gene, exon, meta files, output folder
			 if (String::Compare (celsFileNameTextBox->Text, S"") == 0)
			 {
				 MessageBox::Show (S"Please select a Cel Group File", S"Error");
				 return;
			 }
			 if (String::Compare (geneDataFileTextBox->Text, S"") == 0)
			 {
				 MessageBox::Show (S"Please select a Gene Level Summary File", S"Error");
				 return;
			 }
			 if (String::Compare (probesetDataFileTextBox->Text, S"") == 0)
			 {
				 MessageBox::Show (S"Please select an Exon Level Summary File", S"Error");
				 return;
			 }
			 if (String::Compare (mapFileTextBox->Text, S"") == 0)
			 {
				 MessageBox::Show (S"Please select a Meta Probeset List", S"Error");
				 return;
			 }
			 if (String::Compare (outFolderTextBox->Text, S"") == 0)
			 {
				 MessageBox::Show (S"Please enter an Output Folder", S"Error");
				 return;
			 }
			 // pass string arguments to midas code as char *
			 // mcppfaq_convertstr.cpp
			 // compile with: /clr /LD /EHsc
			 const char* celsFileName =
				 (const char*)(Marshal::StringToHGlobalAnsi(celsFileNameTextBox->Text)).ToPointer();
			 const char* geneDataFileName =
				 (const char*)(Marshal::StringToHGlobalAnsi(geneDataFileTextBox->Text)).ToPointer();
			 const char* probesetDataFileName =
				 (const char*)(Marshal::StringToHGlobalAnsi(probesetDataFileTextBox->Text)).ToPointer();
			 const char* mapFileName =
				 (const char*)(Marshal::StringToHGlobalAnsi(mapFileTextBox->Text)).ToPointer();
			 const char* outputDirectory =
				 (const char*)(Marshal::StringToHGlobalAnsi(outFolderTextBox->Text)).ToPointer();

			 // need to put log stabilization factor conversion inside a try block - might not work
			 bool conversionOk = false;
			 float logStabilize = 8.0;	// default
			 try
			 {
				 logStabilize = System::Convert::ToSingle(logStabilizeTextBox->Text);
				 conversionOk = true;
			 }
			 catch (System::OverflowException*)
			 {
				 String* errMsg = S"The conversion of the log stabilization factor overflowed.";
				 MessageBox::Show (errMsg, S"Error");
			 }
			 catch (System::FormatException*)
			 {
				 String* msg1 = S"The Log Stabilization Factor ";
				 String* msg2 = S" cannot be converted into a float.";
				 String* errMsg = String::Concat (msg1, logStabilizeTextBox->Text, msg2);
				 MessageBox::Show (errMsg, S"Error");
			 }
			 bool outputDirectoryOk = true;
			 if (conversionOk && strlen (outputDirectory) != 0)
			 {
				 // create a reference to a directory
				 DirectoryInfo* di = new DirectoryInfo(outFolderTextBox->Text);
				 // create the directory only if it does not already exist
				 if (di->Exists == false)
				 {
					 try
					 {
						 // make sure di->Create() worked
						 outputDirectoryOk = false;
						 di->Create();
						 outputDirectoryOk = true;
					 }
					 catch (System::ArgumentException*)
					 {
						 String* errMsg = S"The requested output directory does not specify a valid file path.";
						 MessageBox::Show (errMsg, S"Error");
					 }
					 catch (System::IO::IOException*)
					 {
						 String* errMsg = S"The requested output directory cannot be created.";
						 MessageBox::Show (errMsg, S"Error");
					 }
					 catch (System::Security::SecurityException*)
					 {
						 String* msg1 = S"The requested output directory could not be created.";
						 String* msg2 = S"The caller does not have the required permission.";
						 String* errMsg = String::Concat (msg1, msg2);
						 MessageBox::Show (errMsg, S"Error");
					 }
				 }		// end if (di->Exists == false)
			 }			// end if (conversionOk && strlen (outputDirectory) != 0)

			 // set up the "actual" wantPvalues, based on the
			 // state of the object's notWantPvalues
			 bool wantPvalues = true;	// default
			 if (this->notWantPvalues)
				 wantPvalues = false;

			 // generate pseudo command line for output file headers
			 String* noPvaluesString;
			 if (!wantPvalues)
				 noPvaluesString = S" -nopvalues";
			 String* fstatsString;
			 if (this->wantFstats)
				 fstatsString = S" -fstats";
			 String* indicesString;
			 if (this->wantIndices)
				 indicesString = S" -indices";
			 String* noLogString;
			 if (this->noLogTransform)
				 noLogString = S" -no-logtrans";
			 String* commandLineArray[] =
			 {
				 S"midas -cels ", celsFileNameTextBox->Text,
					 S" -genedata ", geneDataFileTextBox->Text,
					 S" -probesetdata ", probesetDataFileTextBox->Text,
					 S" -map ", mapFileTextBox->Text,
					 S" -out-dir ", outFolderTextBox->Text,
					 noPvaluesString, fstatsString, indicesString,
					 S" -stabilize ", logStabilizeTextBox->Text,
					 noLogString
			 };

			 // concatenate string array, convert to C++ string
			 String* commandLineString = String::Concat (commandLineArray);
			 const char* commandLineChars =
				 (const char*)(Marshal::StringToHGlobalAnsi(commandLineString)).ToPointer();
			 std::string commandLine = commandLineChars;

			 if (conversionOk && outputDirectoryOk)
			 {
				 try
				 {
					 // create object to configure, run midas
					 midasConfigureRun configureRun (celsFileName, geneDataFileName,
						 probesetDataFileName, mapFileName, outputDirectory,
						 wantPvalues, this->wantFstats, this->wantIndices,
						 logStabilize, commandLine, this->execVersion, this->noLogTransform);
					 // configure step may return a non-fatal warning message
					 bool okToRun = true;
					 std::string* warningMsgstring = configureRun.configure();
					 if (warningMsgstring != 0)
					 {
						 String* warningMsg (warningMsgstring->c_str());
						 if (MessageBox::Show
							 (warningMsg, S"Warning", MessageBoxButtons::OKCancel)
							    != DialogResult::OK)
						 {
							 // user chose cancel - back out newly created output files
							 configureRun.deleteOutputs();
							 okToRun = false;
						 }
					 }
					 if (okToRun)
					 {
						 // display "processing" box
						 MidasProcessingBox* processing = new MidasProcessingBox();
						 processing->Show();

						 // perform midas run
						 configureRun.run();

						 // run complete
						 // close the processing box, notify user
						 processing->Close();
						 MessageBox::Show (S"Run completed successfully.", S"Finished");
					 }

				 }
				 catch (std::exception &e)
				 {
					 MessageBox::Show (e.what(), S"Error");
				 }
			 }

			 Marshal::FreeHGlobal(System::IntPtr((void*)celsFileName));
			 Marshal::FreeHGlobal(System::IntPtr((void*)geneDataFileName));
			 Marshal::FreeHGlobal(System::IntPtr((void*)probesetDataFileName));
			 Marshal::FreeHGlobal(System::IntPtr((void*)mapFileName));
			 Marshal::FreeHGlobal(System::IntPtr((void*)outputDirectory));
			 Marshal::FreeHGlobal(System::IntPtr((void*)commandLineChars));
		 }

private: System::Void exitButton_Click(System::Object *  sender, System::EventArgs *  e)
		 {
			 Application::Exit();
		 }

private: System::Void label1_Click(System::Object *  sender, System::EventArgs *  e)
		 {
		 }

private: System::Void menuItem2_Click_1(System::Object *  sender, System::EventArgs *  e)
		 {
		 }

private: System::Void indicesGroupBox_Enter(System::Object *  sender, System::EventArgs *  e)
		 {
		 }

private: System::Void pValuesGroupBox_Enter(System::Object *  sender, System::EventArgs *  e)
		 {
		 }

private: System::Void outFolderButton_Click(System::Object *  sender, System::EventArgs *  e)
		 {
			 const char* folderPath =
				 (const char*)(Marshal::StringToHGlobalAnsi(outFolderTextBox->Text)).ToPointer();
			 // open folder dialog box
			 FolderDialog dlg (folderPath);
			 if (dlg.DoModal() == IDOK)
			 {
				 String* pathString (dlg.GetPathName().c_str());
				 outFolderTextBox->Text = pathString;
			 }
			 Marshal::FreeHGlobal (System::IntPtr ((void*)folderPath));
		 }


private: System::Void logTransRadioNo_Click(System::Object *  sender, System::EventArgs *  e)
		 {
			 if (sender == logTransRadioYes)
				 this->noLogTransform = false;
			 else if (sender == logTransRadioNo)
				 this->noLogTransform = true;
		 }

};
}


